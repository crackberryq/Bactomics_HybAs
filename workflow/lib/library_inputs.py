"""Experimental paired-end library adapter; no taxonomy or AWARE policy here.

Only ordinary FR paired-end libraries are executable. Read names are namespaced
before shared classification, then used to restore library boundaries after the
existing sample/platform cleaning decision. Original headers are retained.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import gzip
import html
from itertools import zip_longest
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


ADAPTER = "illumina-libraries-v1"
MODES = {"illumina_only", "dual_illumina_runs", "dual_illumina_complementary"}


def execution_samples(validated):
    """Gate execution separately from the broader metadata-only validator."""
    samples = {}
    for sample, selection in zip(validated["manifest"]["samples"], validated["selections"]):
        name = sample["sample_id"]
        if selection["selected_mode"] not in MODES:
            raise ValueError(f"{name}: manifest mode {selection['selected_mode']!r} is metadata-only; "
                             "this execution adapter supports Illumina paired-end libraries only")
        for lib in sample["libraries"]:
            if lib["library_type"] != "paired_end" or lib["orientation"] != "fr":
                raise ValueError(f"{name}/{lib['library_id']}: only paired_end/fr preprocessing is executable; "
                                 "mate-pair and rf/ff libraries remain metadata-only")
        samples[name] = {**sample, "selection": selection, "adapter": ADAPTER,
                         "assembler": "spades", "qualification": "experimental_software_adapter"}
    return samples


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="ascii") as handle:
        number = 0
        while header := handle.readline():
            number += 1
            seq, plus, qual = [handle.readline().rstrip("\r\n") for _ in range(3)]
            header = header.rstrip("\r\n")
            if not header.startswith("@") or not header[1:].strip() or not plus.startswith("+"):
                raise ValueError(f"{path}: malformed four-line FASTQ record {number}")
            if not seq or len(seq) != len(qual) or any(c.isspace() for c in seq):
                raise ValueError(f"{path}: sequence/quality length mismatch or empty record {number}")
            if any(not 33 <= ord(c) <= 126 for c in qual):
                raise ValueError(f"{path}: invalid quality character at record {number}")
            yield header, seq, qual


def read_id(header, mate):
    fields = header[1:].split()
    name = fields[0]
    if name.endswith(("/1", "/2")):
        if name[-1] != str(mate):
            raise ValueError(f"Wrong mate label in R{mate}: {header}")
        name = name[:-2]
    if len(fields) > 1 and re.match(r"[12]:", fields[1]) and fields[1][0] != str(mate):
        raise ValueError(f"Wrong Illumina mate label in R{mate}: {header}")
    return name


def pairs(r1, r2):
    for number, (left, right) in enumerate(zip_longest(fastq(r1), fastq(r2)), 1):
        if left is None or right is None:
            raise ValueError(f"Unequal paired FASTQ record counts at pair {number}: {r1}, {r2}")
        if read_id(left[0], 1) != read_id(right[0], 2):
            raise ValueError(f"Mismatched paired read IDs at pair {number}: {r1}, {r2}")
        yield left, right


def record(handle, value, name=None):
    header, seq, qual = value
    handle.write(f"{'@' + name if name else header}\n{seq}\n+\n{qual}\n")


def prepare(plan, r1, r2, audit):
    """Fully audit pairs and namespace IDs; raw input files are never modified."""
    audit = Path(audit)
    audit.mkdir(parents=True, exist_ok=True)
    counts = []
    with gzip.open(r1, "wt") as left, gzip.open(r2, "wt") as right, \
            gzip.open(audit / "original_headers.jsonl.gz", "wt") as identities:
        for li, lib in enumerate(plan["libraries"]):
            for gi, group in enumerate(lib["read_groups"]):
                count = 0
                for count, (a, b) in enumerate(pairs(group["r1"], group["r2"]), 1):
                    name = f"hybasL{li}G{gi}R{count}"
                    record(left, a, name)
                    record(right, b, name)
                    identities.write(json.dumps({"id": name, "r1": a[0], "r2": b[0]}) + "\n")
                if not count:
                    raise ValueError(f"{lib['library_id']}/{group['read_group_id']}: no read pairs")
                counts.append({"library_id": lib["library_id"], "read_group_id": group["read_group_id"],
                               "namespace": f"hybasL{li}G{gi}", "pairs": count})
    write_json(audit / "plan.json", plan)
    write_json(audit / "read_groups.json", counts)


def trim(plan, r1, r2, outdir, combined1, combined2, report_html, report_json, threads):
    """Restore libraries after cleaning; run fastp separately for each library."""
    outdir = Path(outdir).absolute()
    outdir.mkdir(parents=True, exist_ok=True)
    counts = [0] * len(plan["libraries"])
    reports, dataset = [], []
    # Clean inputs are scratch files, not permanent second copies of all reads.
    with tempfile.TemporaryDirectory(prefix=".library-clean-", dir=outdir) as temporary:
        paths = [(Path(temporary) / f"{li}_R1.fq.gz", Path(temporary) / f"{li}_R2.fq.gz")
                 for li in range(len(counts))]
        with ExitStack() as stack:
            handles = [(stack.enter_context(gzip.open(a, "wt")), stack.enter_context(gzip.open(b, "wt")))
                       for a, b in paths]
            previous = None
            for a, b in pairs(r1, r2):
                match = re.fullmatch(r"hybasL(\d+)G(\d+)R([1-9]\d*)", read_id(a[0], 1))
                if not match:
                    raise ValueError("Library identity lost during cleaning")
                li, gi, ri = map(int, match.groups())
                key = (li, gi, ri)
                if li >= len(counts) or gi >= len(plan["libraries"][li]["read_groups"]):
                    raise ValueError("Unknown library/read-group identity after cleaning")
                if previous is not None and key <= previous:
                    raise ValueError("Duplicate or reordered read identities after cleaning")
                previous = key
                record(handles[li][0], a)
                record(handles[li][1], b)
                counts[li] += 1
        for li, lib in enumerate(plan["libraries"]):
            name = lib["library_id"]
            if not counts[li]:
                raise ValueError(f"{name}: no pairs remain after cleaning; refusing to silently drop a library")
            left, right = [outdir / f"{name}_R{mate}.trimmed.fq.gz" for mate in (1, 2)]
            html_path, json_path = outdir / f"{name}.fastp.html", outdir / f"{name}.fastp.json"
            subprocess.run(["fastp", "-i", str(paths[li][0]), "-I", str(paths[li][1]),
                            "-o", str(left), "-O", str(right), "--thread", str(threads),
                            "--html", str(html_path), "--json", str(json_path)], check=True)
            remaining = sum(1 for _ in pairs(left, right))
            if not remaining:
                raise ValueError(f"{name}: no pairs remain after fastp; refusing to silently drop a library")
            reports.append({"library_id": name, "library_type": lib["library_type"],
                            "orientation": lib["orientation"], "declared_insert_size_bp": lib.get("insert_size_bp"),
                            "cleaned_pairs": counts[li], "trimmed_pairs": remaining,
                            "fastp": json.loads(json_path.read_text()), "html_path": str(html_path)})
            dataset.append({"type": "paired-end", "orientation": lib["orientation"],
                            "left reads": [str(left)], "right reads": [str(right)]})
    # SPAdes estimates insert distributions. Declared sizes order the libraries,
    # but are not passed as an invented/unsupported dataset attribute.
    order = sorted(range(len(dataset)), key=lambda i: (plan["libraries"][i].get("insert_size_bp", 0), i))
    write_json(outdir / "dataset.yaml", [dataset[i] for i in order])  # JSON is valid YAML.
    write_json(outdir / "libraries.json", {"adapter": ADAPTER, "plan": plan,
               "spades_library_order": [plan["libraries"][i]["library_id"] for i in order], "libraries": reports})
    for mate, dest in (("left reads", combined1), ("right reads", combined2)):
        # Concatenated gzip members retain complete, synchronized paired streams.
        with open(dest, "wb") as handle:
            for lib in dataset:
                with open(lib[mate][0], "rb") as source:
                    shutil.copyfileobj(source, handle)
    write_json(report_json, {"schema": "hybas-library-fastp-v1", "libraries": reports})
    panels = []
    for report in reports:
        panels.append(f"<details><summary>{html.escape(report['library_id'])}: "
                      f"{report['cleaned_pairs']} cleaned pairs, {report['trimmed_pairs']} trimmed pairs</summary>"
                      f"<iframe title=\"fastp report\" style=\"width:100%;height:800px\" "
                      f"srcdoc=\"{html.escape(Path(report['html_path']).read_text(), quote=True)}\"></iframe></details>")
    Path(report_html).write_text("<!doctype html><meta charset=\"utf-8\"><title>Library trimming</title>"
                                "<h1>Illumina library trimming</h1><p>Separate fastp reports; "
                                "rates are not pooled across libraries.</p>" + "".join(panels))


def promote(source, dest, minimum):
    """Promote contigs, applying the same declared minimum-length contract."""
    kept = 0
    header, sequence = None, []
    with open(source) as handle, open(dest, "w") as output:
        def emit():
            nonlocal kept
            seq = "".join(sequence)
            if header and len(seq) >= minimum:
                output.write(header + "\n" + seq + "\n")
                kept += 1
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                emit()
                header, sequence = line, []
            elif line:
                if header is None:
                    raise ValueError("SPAdes output is not FASTA")
                sequence.append(line)
        emit()
    if not kept:
        raise ValueError(f"No SPAdes contigs meet min_fasta_length={minimum}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    prep = sub.add_parser("prepare")
    for arg in ("plan", "r1", "r2", "audit"):
        prep.add_argument("--" + arg, required=True)
    trimming = sub.add_parser("trim")
    for arg in ("plan", "r1", "r2", "outdir", "combined1", "combined2", "report-html", "report-json"):
        trimming.add_argument("--" + arg, required=True)
    trimming.add_argument("--threads", type=int, required=True)
    promotion = sub.add_parser("promote")
    promotion.add_argument("source")
    promotion.add_argument("dest")
    promotion.add_argument("minimum", type=int)
    args = vars(parser.parse_args())
    action = args.pop("action")
    if action == "prepare":
        args["plan"] = json.loads(args["plan"])
    elif action == "trim":
        args["plan"] = json.loads(Path(args["plan"]).read_text())
    {"prepare": prepare, "trim": trim, "promote": promote}[action](**args)


if __name__ == "__main__":
    main()

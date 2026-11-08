"""Pure validation and reporting helpers for the directly runnable HybAs workflow."""

import gzip
import hashlib
import json
import math
import os
import re
import shutil
from pathlib import Path

VERSION = "HybAs-AWARE-10.0.0a1"


def normalize_taxid(value):
    """Normalize numeric IDs without turning YAML null into a manual override."""
    if value is None or value == "":
        return ""
    if isinstance(value, bool):
        raise ValueError("target_taxid must be a numeric taxonomy ID or empty, not a boolean")
    text = str(value).strip()
    if not text:
        return ""
    if not re.fullmatch(r"[0-9]+(?:\.0)?", text):
        raise ValueError(f"Invalid taxonomy ID {value!r}; use an integer or YAML null/empty string")
    return str(int(text.removesuffix(".0")))


def manual_taxid(value):
    result = normalize_taxid(value)
    if result == "0":
        raise ValueError("target_taxid must be positive; 0 denotes unclassified reads")
    return result


def validate_isolate_names(names):
    seen = set()
    for name in names:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) or name in {".", ".."}:
            raise ValueError(f"Invalid isolate identifier {name!r}: use letters, digits, '.', '_' or '-'")
        if name in seen:
            raise ValueError(f"Duplicate isolate identifier: {name}")
        seen.add(name)


def pair_illumina_files(r1, r2):
    """Pair lanes by the filename surrounding the final _R1/_R2 marker."""
    def keyed(paths, mate):
        out = {}
        for path in paths:
            name = Path(path).name
            matches = list(re.finditer(r"_R" + str(mate) + r"(?=[_.])", name))
            if not matches:
                raise ValueError(f"Cannot identify mate {mate} in {name}")
            match = matches[-1]
            key = name[:match.start()] + "_R?" + name[match.end():]
            key = re.sub(r"\.(?:fastq|fq)(?:\.gz)?$", "", key)
            if key in out:
                raise ValueError(f"Duplicate lane/mate input: {out[key]} and {path}")
            out[key] = path
        return out
    left, right = keyed(r1, 1), keyed(r2, 2)
    if left.keys() != right.keys():
        raise ValueError(f"Unpaired Illumina lanes: missing R2={sorted(left.keys()-right.keys())}; "
                         f"missing R1={sorted(right.keys()-left.keys())}")
    return [left[k] for k in sorted(left)], [right[k] for k in sorted(right)]


def validate_input_files(paths):
    for path in paths:
        if not os.path.isfile(path) or not os.access(path, os.R_OK) or os.path.getsize(path) == 0:
            raise ValueError(f"Input is missing, empty or unreadable: {path}")
        opener = gzip.open if str(path).endswith(".gz") else open
        try:
            with opener(path, "rt", encoding="ascii") as handle:
                h, seq, plus, qual = [handle.readline().rstrip("\r\n") for _ in range(4)]
            if not h.startswith("@") or not plus.startswith("+") or not seq or len(seq) != len(qual):
                raise ValueError("invalid first FASTQ record")
        except (OSError, UnicodeError, EOFError, ValueError) as exc:
            raise ValueError(f"Cannot read FASTQ input {path}: {exc}") from exc


def polish_plan(mode, requested, rounds, polypolish):
    stages = []
    if mode in {"hybrid", "ont_only"} and requested != "none":
        stages += [f"racon{i}" for i in range(1, rounds + 1)]
        if requested in {"medaka", "tripolish"}:
            stages.append("medaka")
    if mode in {"hybrid", "illumina_only", "illumina_pacbio"} and requested == "tripolish" and polypolish:
        stages.append("polypolish")
    return stages


def busco_db_ready(path):
    root = Path(path)
    # A directory or an old .complete stamp alone is not a completed dataset.
    return (root.joinpath("dataset.cfg").is_file()
            and root.joinpath("dataset.cfg").stat().st_size > 0
            and any(p.stat().st_size > 0 for p in root.glob("hmms/*.hmm"))
            and root.joinpath("scores_cutoff").is_file()
            and root.joinpath("scores_cutoff").stat().st_size > 0
            and root.joinpath("lengths_cutoff").is_file()
            and root.joinpath("lengths_cutoff").stat().st_size > 0)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path):
    path = Path(path)
    if not path.is_file():
        return {"path": str(path), "exists": False}
    st = path.stat()
    return {"path": str(path.resolve()), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def require_fasta(path):
    count = 0
    length = 0
    current_length = 0
    with open(path, encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if count and not current_length:
                    raise ValueError(f"Empty FASTA record in {path}")
                count += 1
                current_length = 0
            else:
                if not count or not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", line):
                    raise ValueError(f"Invalid assembly FASTA in {path}")
                length += len(line)
                current_length += len(line)
    if not count or not length or not current_length:
        raise ValueError(f"Assembly FASTA is empty: {path}")
    return {"contigs": count, "bases": length, "sha256": sha256_file(path)}


def quast_metric(path, name):
    """QUAST report.tsv has one metric per row (also accept transposed exports)."""
    rows = [line.rstrip("\n").split("\t") for line in Path(path).read_text().splitlines() if line.strip()]
    for row in rows:
        if len(row) > 1 and row[0] == name:
            return row[1]
    if len(rows) > 1 and name in rows[0]:
        index = rows[0].index(name)
        if index < len(rows[1]):
            return rows[1][index]
    return "NA"


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)

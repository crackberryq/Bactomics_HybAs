"""Measure retained polishing states and export the HybAs research trajectory."""
import argparse
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import subprocess
from collections import Counter, defaultdict
from urllib.parse import quote


def write_table(path, rows, fields):
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def fasta_lengths(path, proteins=False):
    lengths = []
    with path.open() as source:
        sequence = []
        for line in source:
            if line.startswith(">"):
                if sequence:
                    lengths.append(length(sequence, proteins))
                sequence = []
            else:
                sequence.append(line.strip())
        if sequence:
            lengths.append(length(sequence, proteins))
    return lengths


def length(parts, proteins):
    seq = "".join(parts)
    return len(seq.removesuffix("*") if proteins else seq)


def markers(path):
    calls = {}
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        marker, status, *_ = line.split("\t")
        if status not in {"Complete", "Duplicated", "Fragmented", "Missing"}:
            raise ValueError(f"Unexpected BUSCO status {status!r}")
        if marker in calls and calls[marker] != status:
            raise ValueError(f"Inconsistent status for {marker}")
        calls[marker] = status
    if not calls:
        raise ValueError(f"Empty BUSCO table: {path}")
    return calls


def components(previous, current):
    return (current["complete_pct"] - previous["complete_pct"],
            previous["short_pct"] - current["short_pct"],
            100 * (current["mean_aa"] - previous["mean_aa"]) / previous["mean_aa"])


def scaled(delta):
    return tuple(x / scale for x, scale in zip(delta, (0.5, 0.1, 0.1)))


def draw(path, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(9, 6), layout="constrained")
    ax = fig.add_subplot(projection="3d")
    x = list(range(len(rows)))
    y = [r["short_pct"] for r in rows]
    z = [r["mean_aa"] for r in rows]
    ax.plot(x, y, z, color="#536e73")
    dots = ax.scatter(x, y, z, c=[r["complete_pct"] for r in rows],
                      cmap="viridis", vmin=0, vmax=100, depthshade=False, s=55)
    ax.set_xticks(x, [r["stage"] for r in rows], rotation=18)
    ax.set_ylabel("Proteins <100 aa (%)")
    ax.set_zlabel("Mean protein length (aa)")
    fig.colorbar(dots, ax=ax, shrink=0.65, label="BUSCO complete (%)")
    fig.savefig(path)
    plt.close(fig)


def feature_index(gff):
    bins = defaultdict(list)
    for line in gff.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        contig, _, kind, start, end, _, _, _, attrs = line.split("\t")
        if kind == "CDS":
            lo, hi = int(start), int(end)
            for block in range(lo // 10000, hi // 10000 + 1):
                bins[contig, block].append((lo, hi, attrs.split(";")[0]))
    return bins


def overlaps(index, contig, pos):
    return ",".join(label for lo, hi, label in index.get((contig, pos // 10000), [])
                    if lo <= pos <= hi)


def edit_evidence(root, pair, before, after):
    ref_features = feature_index(root / before / "genes.gff")
    qry_features = feature_index(root / after / "genes.gff")
    counts = Counter()
    path = root / (pair + ".edits.tsv")
    with path.open("w", newline="") as out, (root / (pair + ".snps")).open() as source:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["ref_contig", "ref_pos", "ref_base", "query_contig", "query_pos",
                         "query_base", "event", "ref_CDS", "query_CDS"])
        for line in source:
            fields = line.split()
            if not fields:
                continue
            rp, rb, qb, qp = fields[:4]
            ref, qry = fields[-2:]
            event = "insertion" if rb == "." else "deletion" if qb == "." else "substitution"
            counts[event] += 1
            writer.writerow([ref, rp, rb, qry, qp, qb, event,
                             overlaps(ref_features, ref, int(rp)) if rb != "." else "",
                             overlaps(qry_features, qry, int(qp)) if qb != "." else ""])
    return counts


def measure(args):
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=True)
    receipt = json.loads(args.receipt.read_text())
    stages = [s for s in receipt["stages"] if s["stage"] not in receipt["skipped_stages"]]
    rows, transitions, changes, links = [], [], [], []
    previous_calls = {}
    previous_asm = None
    previous_q = (0, 0, 0)
    cumulative = peak = 0
    for stage in stages:
        name = stage["stage"]
        asm = Path(stage["path"]).resolve()
        with asm.open("rb") as source:
            digest = hashlib.file_digest(source, "sha256").hexdigest()
        if digest != stage["sha256"]:
            raise ValueError(f"Assembly checksum changed: {name}")
        dest = root / name
        dest.mkdir()
        subprocess.run(["prodigal", "-p", "single", "-g", "11", "-i", asm,
                        "-a", dest / "proteins.faa", "-f", "gff", "-o", dest / "genes.gff"], check=True)
        subprocess.run(["busco", "--offline", "-m", "genome", "-l", args.lineage.resolve(),
                        "-i", asm, "-c", str(args.threads), "--out_path", dest, "-o", "busco"], check=True)
        full, = (dest / "busco").glob("run_*/full_table.tsv")
        calls = markers(full)
        counts = Counter(calls.values())
        proteins = fasta_lengths(dest / "proteins.faa", True)
        lengths = sorted(fasta_lengths(asm), reverse=True)
        total = sum(lengths)
        seen = 0
        for n50 in lengths:
            seen += n50
            if 2 * seen >= total:
                break
        short = sum(n < 100 for n in proteins)
        row = dict(stage=name, sha256=stage["sha256"], fasta=os.path.relpath(asm, root),
                   contigs=len(lengths), bases=total, n50=n50, protein_count=len(proteins),
                   short_count=short, short_pct=100 * short / len(proteins),
                   mean_aa=sum(proteins) / len(proteins))
        for label, count in dict(single=counts["Complete"], duplicated=counts["Duplicated"],
                                 complete=counts["Complete"]+counts["Duplicated"],
                                 fragmented=counts["Fragmented"], missing=counts["Missing"]).items():
            row[label] = count
            row[label+"_pct"] = 100 * count / len(calls)
        links.extend([(name+" assembly", row["fasta"]), (name+" BUSCO markers", str(full.relative_to(root))),
                      (name+" proteins", name+"/proteins.faa"), (name+" features", name+"/genes.gff")])
        if rows:
            previous = rows[-1]
            delta = components(previous, row)
            z = scaled(delta)
            votes = [1 if x >= 1-1e-10 else -1 if x <= -1+1e-10 else 0 for x in z]
            label = "R" if votes.count(1)>=2 else "D" if votes.count(-1)>=2 else "M"
            pair = previous["stage"]+"-to-"+name
            q = scaled(components(rows[0], row))
            cumulative += math.dist(previous_q, q)
            previous_q = q
            peak = max(peak, math.hypot(*q))
            if calls.keys() != previous_calls.keys():
                raise ValueError("BUSCO marker sets differ between stages")
            for marker in sorted(calls):
                if calls[marker] != previous_calls[marker]:
                    changes.append(dict(pair=pair, marker=marker, before=previous_calls[marker], after=calls[marker]))
            subprocess.run(["dnadiff", "-p", root / pair, previous_asm, asm], check=True)
            edits = edit_evidence(root, pair, previous["stage"], name)
            transitions.append(dict(pair=pair, delta_busco_pp=delta[0], delta_short_pp=-delta[1],
                                    relative_mean_change_pct=delta[2], z_busco=z[0], z_short=z[1],
                                    z_length=z[2], votes=",".join(map(str,votes)), response=label,
                                    aligned_substitutions=edits["substitution"], aligned_inserted_bp=edits["insertion"],
                                    aligned_deleted_bp=edits["deletion"]))
            links.extend([(pair+" edit locations and overlapping CDS", pair+".edits.tsv"), (pair+" alignment and unaligned bases", pair+".report")])
        rows.append(row)
        previous_calls, previous_asm = calls, asm
    write_table(root/"stages.tsv", rows, rows[0].keys())
    if transitions:
        write_table(root/"transitions.tsv", transitions, transitions[0].keys())
    write_table(root/"marker_changes.tsv", changes, ["pair","marker","before","after"])
    (root/"geometry.json").write_text(json.dumps(dict(cumulative_path=cumulative, endpoint_displacement=math.hypot(*previous_q),
                                                     peak_displacement=peak, baseline=rows[0]["stage"]), indent=2))
    versions = {}
    for tool, flag in [("busco","--version"), ("prodigal","-v"), ("nucmer","--version")]:
        versions[tool] = subprocess.check_output([tool, flag], stderr=subprocess.STDOUT, text=True).strip()
    (root/"recipe.json").write_text(json.dumps(dict(versions=versions,lineage=str(args.lineage.resolve()),
        dataset=(args.lineage/"dataset.cfg").read_text(), prodigal_mode="single", genetic_code=11,
        remove_one_terminal_stop=True, thresholds=[0.5,0.1,0.1], matched_control=None), indent=2))
    mapping = args.receipt.resolve().parents[1] / "work/mapping"
    for filename in ("assembly.medaka.fa", "illumina.bam", "illumina.bam.bai",
                     "polypolish_R1.sam", "polypolish_R2.sam"):
        evidence = mapping / filename
        if evidence.is_file():
            links.append(("Polypolish input mapping: " + filename, os.path.relpath(evidence, root)))
    draw(root/"trajectory.svg", rows)
    def display(items, fields):
        headers = "".join(f"<th>{html.escape(key)}</th>" for key in fields)
        body = []
        for item in items:
            values = [format(item[key], ".6g") if isinstance(item[key], float) else str(item[key]) for key in fields]
            body.append("<tr>" + "".join(f"<td>{html.escape(value)}</td>" for value in values) + "</tr>")
        return "<table><tr>" + headers + "</tr>" + "".join(body) + "</table>"

    summary = display(rows, ["stage", "complete_pct", "short_count", "protein_count", "short_pct", "mean_aa", "n50"])
    if transitions:
        summary += display(transitions, ["pair", "votes", "response", "aligned_substitutions", "aligned_inserted_bp", "aligned_deleted_bp"])
    summary += (f"<p>Cumulative path {cumulative:.3f} · endpoint displacement {math.hypot(*previous_q):.3f} · "
                f"peak displacement {peak:.3f} (threshold units relative to the initial assembly).</p>")
    index = [(p.name,p.name) for p in sorted(root.iterdir()) if p.suffix in {".tsv",".json",".svg"}]
    (root/"index.html").write_text('<meta charset="utf-8"><title>Polishing trajectory</title><h1>Polishing trajectory</h1>'
        '<style>body{font:16px system-ui;max-width:1100px;margin:2rem auto;color:#16363b}td,th{padding:.4rem;text-align:left;border-bottom:1px solid #ccc}img{max-width:100%}</style>'
        '<img alt="Ordered polishing trajectory" width="800" src="trajectory.svg"><p>Stage spacing represents order. '
        'Transitions describe gene-space change. Edit counts cover one-to-one aligned sequence. '
        'Unaligned sequence is summarized in each alignment report.</p>'+summary+'<ul>'+
        ''.join(f'<li><a href="{quote(url)}">{html.escape(title)}</a></li>' for title,url in index+links)+'</ul>')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--lineage", required=True, type=Path)
    parser.add_argument("--threads", required=True, type=int)
    measure(parser.parse_args())

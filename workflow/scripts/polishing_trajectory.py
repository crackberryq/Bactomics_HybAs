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
import re
import sys
import base64
from collections import Counter, defaultdict
from urllib.parse import quote

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from report_theme import evidence_page
from polishing_validation_figures import projections


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
    rows, transitions, changes = [], [], []
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
                        "-i", asm, "-c", str(args.threads), "--out_path", dest, "-o", "busco"], check=True, cwd=dest)
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
    draw(root / "trajectory.svg", rows)
    render_saved(root, root, args.validation_dir)


def read_table(path):
    with path.open(encoding="utf8", newline="") as source:
        return list(csv.DictReader(source, delimiter="\t"))


def metric_figure(rows, key, title, unit, color):
    """Use a common stage axis and an independent labelled scale for each measure."""
    values = [float(row[key]) for row in rows]
    width = max(445, len(rows) * 75 + 50)
    low, high = min(values), max(values)
    points, marks = [], []
    for i, (row, value) in enumerate(zip(rows, values)):
        fraction = i / (len(rows) - 1) if len(rows) > 1 else 0.5
        x = 45 + (width - 90) * fraction
        y = 101 - 66 * ((value - low) / (high - low) if high != low else 0.5)
        points.append(f"{x:.1f},{y:.1f}")
        label = f"{value:,.0f}" if key == "protein_count" else f"{value:,.2f}"
        marks.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.3" fill="{color}"/>'
                     f'<text x="{x:.1f}" y="{y-10:.1f}">{label}</text>'
                     f'<text x="{x:.1f}" y="133">{html.escape(row["stage"])}</text>')
    return (f'<figure class="evidence-chart"><figcaption><strong>{title}</strong> ({unit})</figcaption>'
            f'<svg viewBox="0 0 {width} 150" role="img" aria-label="{title} by stage">'
            f'<path d="M35 112H{width-35}" stroke="#dce4e3"/><polyline fill="none" stroke="{color}" '
            f'stroke-width="2.4" points="{" ".join(points)}"/><g text-anchor="middle" font-size="11" '
            f'fill="#243f51">{"".join(marks)}</g></svg></figure>')


def evidence_table(items, fields):
    heads = "".join(f"<th scope='col'>{label}</th>" for key, label in fields)
    cells = []
    for item in items:
        values = []
        for key, _ in fields:
            raw = item.get(key, "")
            if raw in (None, ""):
                value = "Unavailable"
            elif key in {"mean_aa", "complete_pct", "short_pct"}:
                value = f"{float(raw):,.2f}"
            elif key == "response":
                value = {"R": "Recovery", "D": "Decline", "M": "Mixed or stable"}.get(raw, raw)
            elif key == "pair":
                value = raw.replace("-to-", " → ")
            else:
                value = f"{int(raw):,}" if str(raw).isdigit() else str(raw)
            values.append(f"<td>{html.escape(value)}</td>")
        cells.append("<tr>" + "".join(values) + "</tr>")
    return '<div class="table-wrap"><table><thead><tr>'+heads+'</tr></thead><tbody>'+''.join(cells)+'</tbody></table></div>'


def render_saved(source, out, validation_dir=None):
    source, out = Path(source).resolve(), Path(out).resolve()
    rows = read_table(source / "stages.tsv")
    if not rows or any(not re.fullmatch(r"[A-Za-z0-9_-]+", row["stage"]) for row in rows):
        raise ValueError("No valid retained stages in stages.tsv")
    pairs = read_table(source / "transitions.tsv") if len(rows) > 1 else []
    keys = [a["stage"]+"-to-"+b["stage"] for a,b in zip(rows, rows[1:])]
    if [r["pair"] for r in pairs] != keys:
        raise ValueError("Transition rows do not match retained stage order")
    out.mkdir(parents=True, exist_ok=True)

    def links(items):
        links_html = []
        for path, label in items:
            if path.is_file():
                href = quote(os.path.relpath(path, out))
                links_html.append(f'<li><a href="{href}">{html.escape(label)}</a></li>')
        return '<ul class="evidence-links">' + "".join(links_html) + '</ul>'

    first, last = rows[0], rows[-1]
    lowest = min(rows, key=lambda row: float(row["complete_pct"]))
    stats = [("Initial BUSCO", f'{float(first["complete_pct"]):.2f}%'),
             ("Lowest BUSCO", f'{float(lowest["complete_pct"]):.2f}%'),
             ("Final retained BUSCO", f'{float(last["complete_pct"]):.2f}%'),
             ("Retained stages", str(len(rows)))]
    body = '<section id="stage"><h2>Effect across stages</h2><div class="grid">' + ''.join(
        f'<div class="stat"><small>{k}</small><strong>{v}</strong></div>' for k,v in stats) + '</div>'
    delta = float(last["complete_pct"]) - float(first["complete_pct"])
    dip = ""
    if lowest is not first and lowest is not last:
        dip = f'The lowest stage was {html.escape(lowest["stage"])} at {float(lowest["complete_pct"]):.2f}%. '
    body += (f'<p class="panel">Complete BUSCO markers changed by <strong>{delta:+.2f} percentage points</strong> '
             f'from {html.escape(first["stage"])} to {html.escape(last["stage"])}. {dip}'
             'Protein predictions and marker recovery describe gene-space effects. Base accuracy requires independent evidence.</p>')
    body += '<div class="chart-grid">' + ''.join(metric_figure(rows,*metric) for metric in [
        ("complete_pct","Complete BUSCO markers","%","#087c78"), ("protein_count","Predicted proteins","count","#2e668e"),
        ("short_pct","Proteins under 100 aa","%","#855e9c"), ("mean_aa","Mean protein length","aa","#986139")])+'</div>'
    body += '<details><summary>Stage measurements</summary>' + evidence_table(rows,[
        ("stage","Stage"),("bases","Assembly (bp)"),("contigs","Contigs"),("n50","N50 (bp)"),("complete_pct","Complete BUSCO (%)"),
        ("fragmented","Fragmented BUSCO"),("missing","Missing BUSCO"),("protein_count","Proteins"),("short_count","Under 100 aa"),("mean_aa","Mean length (aa)")]) + '</details></section>'
    body += '<section id="transitions"><h2>Adjacent-stage changes</h2><div class="panel">'
    if pairs:
        body += ('<p>Each row compares consecutive assemblies across their one-to-one aligned sequence. '
                 'Insertions and deletions are counted in bases.</p>' + evidence_table(pairs,[
                     ("pair","Transition"),("response","Gene-space response"),("aligned_substitutions","Substitutions"),
                     ("aligned_inserted_bp","Inserted (bp)"),("aligned_deleted_bp","Deleted (bp)")]))
        body += ('<p class="note">The response follows the majority direction across BUSCO completeness, the fraction of short proteins and mean protein length. The recipe records the thresholds.</p>')
    else:
        body += '<p>A single retained stage has no adjacent-stage comparison.</p>'
    validation_page = (Path(validation_dir) if validation_dir else source.parent / "polishing_validation") / "index.html"
    if validation_dir or validation_page.is_file():
        href = quote(os.path.relpath(validation_page, out))
        body += f'<p><a href="{href}#edit-density-atlas">See where edits occurred</a> in final assembly coordinates.</p>'
    body += '</div></section><section id="geometry"><h2>Geometric trajectory</h2><div class="panel">'
    body += '<h3>Top, front and side views</h3>' + projections(rows)
    svg = source / "trajectory.svg"
    if svg.is_file():
        image = base64.b64encode(svg.read_bytes()).decode("ascii")
        body += ('<details><summary>Oblique view</summary><img alt="Ordered polishing trajectory" '
                 f'src="data:image/svg+xml;base64,{image}"><p>Stage spacing shows order. '
                 'Colour represents BUSCO completeness.</p></details>')
    if (source / "geometry.json").is_file():
        geo = json.loads((source / "geometry.json").read_text())
        body += (f'<p>In threshold units, the cumulative path is {geo["cumulative_path"]:.3f}, '
                 f'the greatest distance from the initial stage is {geo["peak_displacement"]:.3f} '
                 f'and the endpoint distance is {geo["endpoint_displacement"]:.3f}.</p>')
    body += '</div></section><section id="files"><h2>Evidence files</h2><div class="panel">'
    body += links((source/f,label) for f,label in [
        ("stages.tsv","Stage measurements"),("transitions.tsv","Transition measurements"),("marker_changes.tsv","Marker changes"),
        ("geometry.json","Trajectory geometry"),("trajectory.svg","Trajectory SVG"),("recipe.json","Methods and versions")])
    for row in rows:
        stage = row["stage"]
        items = [(source/row["fasta"],"Assembly FASTA"),(source/stage/"proteins.faa","Protein FASTA"),(source/stage/"genes.gff","Predicted features")]
        items.extend((p,"BUSCO marker table") for p in sorted((source/stage/"busco").glob("run_*/full_table.tsv")))
        body += f'<details><summary>{html.escape(stage)} evidence</summary>{links(items)}</details>'
    for pair in keys:
        body += f'<details><summary>{html.escape(pair.replace("-to-"," → "))} alignment</summary>' + links([
            (source/(pair+".edits.tsv"),"Edits and overlapping CDS"),(source/(pair+".report"),"Alignment and unaligned bases")]) + '</details>'
    mapping = source.parent.parent / "work/mapping"
    files = [(mapping/name,name) for name in ("assembly.medaka.fa","illumina.bam","illumina.bam.bai","polypolish_R1.sam","polypolish_R2.sam") if (mapping/name).is_file()]
    if files:
        body += '<details><summary>Polypolish input mapping</summary>'+links(files)+'</details>'
    body += '</div></section>'
    page = evidence_page("Polishing trajectory", " → ".join(r["stage"] for r in rows), body,
                         [("stage","Stage effects"),("transitions","Adjacent changes"),("geometry","Trajectory views"),("files","Evidence files")],
                         report_type="trajectory")
    (out / "index.html").write_text(page, encoding="utf8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--render-from", type=Path, help="Render retained stage and transition tables")
    parser.add_argument("--validation-dir", type=Path, help="Validation report generated alongside this trajectory")
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--lineage", type=Path)
    parser.add_argument("--threads", type=int)
    args = parser.parse_args()
    if args.render_from:
        render_saved(args.render_from, args.out, args.validation_dir)
    elif args.receipt and args.lineage and args.threads and args.threads > 0:
        measure(args)
    else:
        parser.error("Measuring requires --receipt, --lineage and a positive --threads")

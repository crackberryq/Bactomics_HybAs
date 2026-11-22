"""Backfill final-coordinate evidence from retained stages without re-running polishing."""
import argparse
import base64
import re
import csv
import fcntl
import gzip
import hashlib
import html
import io
import json
import os
import shutil
import subprocess
import tempfile
from bisect import bisect_left
from collections import Counter, defaultdict
from pathlib import Path

from polishing_trajectory import fasta_lengths, write_table
from polishing_validation_figures import density_atlas, projections, spatial_context, complexity_plot
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'lib'))
from report_theme import evidence_page
from urllib.parse import quote

WINDOW = 10000


def records(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def sequences(lines):
    seqs = {}
    name = None
    for line in lines:
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in seqs:
                raise ValueError(f"Duplicate contig {name}")
            seqs[name] = []
        elif line.strip():
            if name is None:
                raise ValueError("Missing FASTA header")
            seqs[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in seqs.items()}


def projection(delta, requested):
    """Map only requested query loci through one-to-one NUCMER blocks."""
    mapped = {}
    with delta.open() as stream:
        next(stream)
        next(stream)
        for line in stream:
            if line.startswith(">"):
                ref, query = line[1:].split()[:2]
                continue
            r, rend, q, qend, *_ = map(int, line.split())
            rd = 1 if rend >= r else -1
            qd = 1 if qend >= q else -1
            coords = requested.get(query, [])
            for gap in stream:
                gap = int(gap)
                size = abs(gap)-1 if gap else (qend-q)*qd+1
                end = q + qd * size
                low = q if qd == 1 else end+1
                high = end if qd == 1 else q+1
                for pos in coords[bisect_left(coords,low):bisect_left(coords,high)]:
                    key = query, pos
                    hit = ref, r + (pos-q)*qd*rd
                    mapped[key] = hit if key not in mapped or mapped[key] == hit else None
                r += rd*size
                q = end
                if not gap:
                    break
                if gap > 0:
                    r += rd
                else:
                    q += qd
    return mapped


def annotation(path, reference):
    counts, features = Counter(), defaultdict(list)
    if not path.is_file() or not path.stat().st_size:
        return counts, features, "Final GFF unavailable"
    intervals = []
    with path.open() as stream:
        for line in stream:
            if line.strip() == "##FASTA":
                break
            cols = line.rstrip().split("\t")
            if len(cols) == 9 and not line.startswith("#"):
                counts[cols[2]] += 1
                if cols[2] == "rRNA":
                    intervals.append((cols[0], int(cols[3]), int(cols[4])))
        embedded = sequences(stream)
    names = defaultdict(list)
    for name, seq in reference.items():
        names[hashlib.sha256(seq.encode()).digest()].append(name)
    aliases = {name: names[hashlib.sha256(seq.encode()).digest()] for name,seq in embedded.items()}
    if not embedded or any(len(x)!=1 for x in aliases.values()) or sorted(x[0] for x in aliases.values()) != sorted(reference):
        return counts, features, "GFF sequence identity unresolved"
    for name, start, end in intervals:
        if name not in aliases or not 1 <= start <= end <= len(reference[aliases[name][0]]):
            raise ValueError("Invalid rRNA coordinates")
        features[aliases[name][0]].append((start,end))
    return counts, features, "Final GFF verified against final FASTA"


def overlap_bp(intervals, lo, hi):
    covered, last = 0, lo-1
    for start,end in sorted(intervals):
        left, right = max(start,lo,last+1), min(end,hi)
        if left <= right:
            covered += right-left+1
            last = right
    return covered


def site_clusters(sites):
    """Bound seek/decompression overhead for nearby loci while limiting memory."""
    positions = []
    contig = None
    for name,pos in sorted(sites):
        if positions and (name!=contig or pos-positions[-1]>500 or pos-positions[0]>=32768):
            yield contig,positions
            positions = []
        contig = name
        positions.append(pos)
    if positions:
        yield contig,positions


def read_support(root, final, iso, sites, threads):
    reads = [iso/"trimmed"/f"R{i}.trimmed.fq.gz" for i in (1,2)]
    fields = ["contig","position","reference_base","informative_bases","A","C","G","T","deletion_spans","anchored_indel_counts"]
    rows = []
    if not sites:
        status = "No projected sites"
    elif not all(p.is_file() and p.stat().st_size for p in reads):
        status = "Paired Illumina reads unavailable"
    else:
        for path in reads:
            with gzip.open(path,"rb") as stream:
                if not stream.read(1):
                    write_table(root/"site_support.tsv", [], fields)
                    return "Paired Illumina reads empty"
        import pysam
        indexed = root/"read_reference.fa"
        shutil.copyfile(final,indexed)
        subprocess.run(["bwa","index",str(indexed)],check=True)
        bam = root/"illumina.final.bam"
        with subprocess.Popen(["bwa","mem","-t",str(threads),str(indexed),*map(str,reads)],stdout=subprocess.PIPE) as bwa:
            try:
                subprocess.run(["samtools","sort","-o",str(bam),"-"],stdin=bwa.stdout,check=True)
            finally:
                bwa.stdout.close()
            if bwa.wait():
                raise RuntimeError("BWA alignment failed")
        pysam.index(str(bam))
        pysam.faidx(str(indexed))
        with pysam.FastaFile(str(indexed)) as fa, pysam.AlignmentFile(str(bam)) as align:
            for contig, positions in site_clusters(sites):
                wanted = set(positions)
                observed = {}
                for col in align.pileup(contig,positions[0]-1,positions[-1],truncate=True,stepper="samtools",min_mapping_quality=20,min_base_quality=20,max_depth=1000000):
                    pos = col.reference_pos + 1
                    if pos not in wanted:
                        continue
                    counts, indels = Counter(), Counter()
                    observed[pos] = counts, indels
                    for hit in col.pileups:
                        if hit.is_refskip:
                            continue
                        if hit.is_del:
                            counts["del"] += 1
                            continue
                        idx = hit.query_position
                        seq = hit.alignment.query_sequence.upper()
                        counts[seq[idx]] += 1
                        if hit.indel > 0:
                            indels["+"+seq[idx+1:idx+1+hit.indel]] += 1
                        elif hit.indel < 0:
                            indels["-"+fa.fetch(contig,pos,pos-hit.indel).upper()] += 1
                for pos in positions:
                    counts, indels = observed.get(pos,(Counter(),Counter()))
                    vals = [counts[b] for b in "ACGT"]
                    rows.append(dict(zip(fields,[contig,pos,fa.fetch(contig,pos-1,pos).upper(),sum(vals),*vals,counts["del"],json.dumps(indels,sort_keys=True)])))
        status = "Final-reference Illumina pileup, mapping and base quality ≥20"
    write_table(root/"site_support.tsv",rows,fields)
    return status


def run(iso, root, threads):
    trajectory = iso/"reports/polishing_trajectory"
    stages = records(trajectory/"stages.tsv")
    if not stages or any(not re.fullmatch(r"[A-Za-z0-9_-]+",s["stage"]) for s in stages):
        raise ValueError("Missing or invalid retained stage identity")
    final = iso/"work/final/assembly.final.fasta"
    with final.open() as stream:
        reference = sequences(stream)
    edits, needed = [], defaultdict(lambda: defaultdict(set))
    for before, after in zip(stages,stages[1:]):
        pair = before["stage"]+"-to-"+after["stage"]
        for row in records(trajectory/(pair+".edits.tsv")):
            side = "query" if row["event"]=="insertion" else "ref"
            stage = after["stage"] if side=="query" else before["stage"]
            contig, pos = row[side+"_contig"], int(row[side+"_pos"])
            needed[stage][contig].add(pos)
            edits.append((pair,stage,contig,pos,row["event"]))
    maps, direct, bins, sites = {}, [], Counter(), set()
    for row in stages:
        stage = row["stage"]
        asm = (trajectory/row["fasta"]).resolve()
        if hashlib.sha256(asm.read_bytes()).hexdigest()!=row["sha256"]:
            raise ValueError(f"Stage checksum changed: {stage}")
        prefix = root/("final-vs-"+stage)
        subprocess.run(["dnadiff","-p",str(prefix),str(final),str(asm)],check=True)
        coords = {k:sorted(v) for k,v in needed[stage].items()}
        maps[stage] = projection(prefix.with_suffix(".1delta"),coords)
        counts = Counter()
        with prefix.with_suffix(".snps").open() as stream:
            for line in stream:
                _,ref,query,*_ = line.split()
                counts["stage_only_bp" if ref=="." else "final_only_bp" if query=="." else "substitutions"] += 1
        lengths = fasta_lengths(trajectory/stage/"proteins.faa",True)
        if len(lengths)!=int(row["protein_count"]) or sum(n<100 for n in lengths)!=int(row["short_count"]):
            raise ValueError(f"Protein counts inconsistent: {stage}")
        direct.append(dict(stage=stage,complete_pct=row["complete_pct"],proteins=len(lengths),under_100=int(row["short_count"]),under_200=sum(n<200 for n in lengths),**{k:counts[k] for k in ("substitutions","stage_only_bp","final_only_bp")}))
    write_table(root/"stage_to_final.tsv",direct,direct[0].keys())
    fields = ["transition","stage","source_contig","source_pos","event","final_contig","final_pos","status"]
    mapped_rows = []
    for pair,stage,contig,pos,event in edits:
        mapped = maps[stage].get((contig,pos))
        if mapped:
            sites.add(mapped)
            bins[mapped[0],(mapped[1]-1)//WINDOW] += 1
        mapped_rows.append(dict(zip(fields,[pair,stage,contig,pos,event,*(mapped or ("","")),"mapped" if mapped else "unresolved"])))
    write_table(root/"projected_edits.tsv",mapped_rows,fields)
    counts, rna, gff_status = annotation(iso/"annotation"/(iso.name+".gff"),reference)
    windows = []
    for contig,seq in reference.items():
        for offset in range(0,len(seq),WINDOW):
            end = min(offset+WINDOW,len(seq))
            windows.append(dict(contig=contig,start=offset+1,end=end,edit_bases=bins[contig,offset//WINDOW],rRNA_overlap_bp=overlap_bp(rna[contig],offset+1,end) if gff_status.startswith("Final GFF verified") else ""))
    write_table(root/"windows.tsv",windows,windows[0].keys())
    hotspots = sorted((w for w in windows if w["edit_bases"]),key=lambda w:(-w["edit_bases"],w["contig"],w["start"]))[:20]
    write_table(root/"hotspots.tsv",hotspots,windows[0].keys())
    # Also observe the preceding anchor for indels represented by pileup after a base.
    evidence_sites = sites | {(c,p-1) for c,p in sites if p>1}
    support = read_support(root,final,iso,evidence_sites,threads)
    meta = dict(final_sha256=hashlib.sha256(final.read_bytes()).hexdigest(),stages={s["stage"]:s["sha256"] for s in stages},
                projected=sum(bins.values()),unresolved=len(edits)-sum(bins.values()),window_bp=WINDOW,
                hotspot_rule="Top 20 nonempty windows by affected base count; ties by contig then start; 1-based inclusive",
                annotation=gff_status,annotation_counts=counts,site_support=support)
    (root/"validation.json").write_text(json.dumps(meta,indent=2))
    render(iso, root, root)


def table(rows):
    if not rows:
        return "<p>No observations.</p>"
    headers = "".join(f"<th>{html.escape(k.replace('_',' '))}</th>" for k in rows[0])
    values = "".join("<tr>" + "".join(f'<td>{html.escape(str(v)) if v != "" else "—"}</td>' for v in row.values()) + "</tr>" for row in rows)
    return f'<div class="table-wrap"><table><thead><tr>{headers}</tr></thead><tbody>{values}</tbody></table></div>'


def atomic_write(path, content):
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False, encoding="utf-8") as out:
        out.write(content)
    os.replace(out.name,path)


def refresh_section(original, fragment, target):
    """Replace the exact validation layouts emitted by the a9 and current writers."""
    if "</main>" not in original:
        raise ValueError("Missing research report main section")
    start = "<!-- HybAs final-coordinate validation -->"
    end = "<!-- End HybAs final-coordinate validation -->"
    section = (start + '<section id="polishing-validation">' + fragment +
               f'<p><a href="{html.escape(target,quote=True)}">Open polishing validation and evidence tables</a></p></section>' + end)
    inserted = False

    def replace(match):
        nonlocal inserted
        if inserted:
            return ""
        inserted = True
        return section

    # Use one pass: a9 could contain both its original inline summary and the backfill.
    old_inline = (r'(?:<section id="polish-overview">)?<h2>Polishing validation</h2>'
                  r'(?:(?!</main>|<!-- HybAs final-coordinate validation -->).)*?'
                  r'<p><a href="\.\./polishing_validation/index\.html">Download polishing validation tables</a></p>')
    marked = (re.escape(start) + r'(?:(?!' + re.escape(start) + r').)*?(?:'
              r'</section>\s*<!-- [Ee]nd HybAs final-coordinate validation -->|'
              r'>Validation evidence tables</a></section>(?:\s*<!-- [Ee]nd HybAs final-coordinate validation -->)?)')
    link_only = (r'<section id="polishing-validation"><h2>Polishing validation</h2><p>'
                 r'<a href="[^"]+">Open polishing validation</a></p></section>')
    pattern = re.compile(marked + "|" + old_inline + "|" + link_only, re.DOTALL)
    updated = pattern.sub(replace, original)
    # Detect a broken marked section before touching the report on disk.
    if updated.count(start) != int(inserted):
        raise ValueError("Incomplete polishing validation section")
    return updated if inserted else updated.replace("</main>",section+"</main>",1)


def render(iso, tables, dest):
    """Render measured tables without invoking native tools or reading raw libraries."""
    traj = iso / "reports/polishing_trajectory"
    stages = records(traj/"stages.tsv")
    meta = json.loads((tables/"validation.json").read_text())
    direct = records(tables/"stage_to_final.tsv")
    names = [row["stage"] for row in stages]
    if not names or names != [row["stage"] for row in direct] or meta["stages"] != {r["stage"]:r["sha256"] for r in stages}:
        raise ValueError("Validation and trajectory stage identities differ")
    final_data = (iso/"work/final/assembly.final.fasta").read_bytes()
    if hashlib.sha256(final_data).hexdigest()!=meta["final_sha256"]:
        raise ValueError("Validation final FASTA checksum differs from retained final assembly")
    reference = sequences(io.StringIO(final_data.decode("utf-8")))
    windows = records(tables/"windows.tsv")
    edits = records(tables/"projected_edits.tsv")
    if sum(int(w["edit_bases"]) for w in windows) != meta["projected"] or sum(r["status"]=="mapped" for r in edits) != meta["projected"] or len(edits)!=meta["projected"]+meta["unresolved"]:
        raise ValueError("Inconsistent projected edit and window totals")
    hotspots = records(tables/"hotspots.tsv")
    observed = {(r["final_contig"],int(r["final_pos"])) for r in edits if r["status"]=="mapped"}
    support = records(tables/"site_support.tsv")
    pairs = [before+"-to-"+after for before,after in zip(names,names[1:])]
    complexity, dispersion, motifs, context = spatial_context(reference, windows, hotspots, edits, support, pairs, int(meta["window_bp"]))
    counts = [int(r["informative_bases"]) for r in support if (r["contig"],int(r["position"])) in observed]
    summary = [
        ("Complete BUSCO at final measured stage",f'{float(stages[-1]["complete_pct"]):.2f}%'),
        ("Projected edited bases",f'{meta["projected"]:,}'),
        ("Excluded from spatial plot",f'{meta["unresolved"]:,}'),
        ("Positions with A/C/G/T observations",f'{sum(v>0 for v in counts):,} of {len(counts):,}' if counts else "Unavailable"),
    ]
    cards = "".join(f'<div class="stat"><small>{label}</small><strong>{value}</strong></div>' for label,value in summary)
    body = '<section id="polish-overview"><h2>Polishing validation</h2><div class="grid">'+cards+'</div></section>'
    body += ('<section id="edit-density-atlas"><h3>Where did polishing change the assembly?</h3><div class="panel">'
             f'<p>The top track counts projected edited bases from all transitions in {int(meta["window_bp"]):,} bp windows. Teal blocks mark windows overlapping rRNA. '
             'Follow a peak down to the transition rows to identify the contributing polishing steps. Repeated edits across steps contribute repeatedly to window totals.</p>'
             + density_atlas(windows, edits, stages) + '</div></section>')
    body += ('<section id="polish-context"><h3>How are edits distributed across sequence contexts?</h3><div class="panel">'
             '<h4>Window complexity and edit burden</h4>' + complexity_plot(complexity) +
             '<p>Each point represents a full window, including windows with zero edits. Entropy measures diversity among overlapping A/C/G/T triplets on the recorded forward strand. '
             'Windows without valid triplets have no entropy estimate. The spatial atlas and downloadable window table retain shorter terminal windows.</p>'
             '<details><summary>Spatial dispersion by polishing step</summary>' + table(dispersion) +
             '<p>The sample variance divided by the mean describes variation in edit counts among full windows. A value above one indicates greater variation than an equal-mean Poisson model. '
             'The ratio requires at least two full windows and a positive mean.</p></details>'
             '<details><summary>Hotspot coordinates, RNA overlap and read depth</summary>' + table(context) +
             '<p>Hotspots are ranked by projected edited bases. RNA overlap counts the union of annotated rRNA bases within each window. '
             'Depth summarizes passing A/C/G/T observations at distinct projected edit positions. A dash denotes unavailable depth. Repeat-family classifications require independent repeat annotations.</p></details>')
    best = {}
    for motif in motifs:
        rank = motif["rank"]
        if rank not in best or motif["frequency_ratio"] > best[rank]["frequency_ratio"]:
            best[rank] = motif
    body += ('<details><summary>Sequence motifs in hotspots</summary>'
             + (table([{k:r[k] for k in ("rank","motif","window_count","background_count","frequency_ratio")} for r in best.values()]) if best else '<p>No enriched triplets with measurable background.</p>')
             + '<p>The strongest frequency ratio for each ranked window is shown. The denominator pools all windows outside the ranked hotspots, including partial windows. '
             'Counts exclude ambiguous triplets and apply a half-count per possible triplet to the frequency estimate. Ratios describe composition. '
             'The export includes every observed triplet with a ratio above one.</p></details></div></section>')
    body += ('<section id="polish-outcomes"><h3>How did each stage differ from the final sequence?</h3><div class="panel">'
             + table(direct) + '<p>Stage-only and final-only columns count bases in one-to-one alignments. Substitutions count aligned mismatches.</p>'
             '<details><summary>Polishing trajectory projections</summary>'+projections(stages))
    oblique = traj/"trajectory.svg"
    if oblique.is_file() and oblique.stat().st_size < 2_000_000:
        data = base64.b64encode(oblique.read_bytes()).decode()
        body += f'<img style="max-width:720px;width:100%" alt="Oblique trajectory coloured by complete BUSCO percentage" src="data:image/svg+xml;base64,{data}">'
    body += '</details></div></section>'
    annotations = meta.get("annotation_counts",{})
    body += '<section id="polish-support"><h3>What evidence supports these changes?</h3><div class="panel"><p>'+html.escape(meta["site_support"])+'.</p>'
    if counts:
        body += (f'<p>{sum(v>0 for v in counts):,} of {len(counts):,} measured edit positions have at least one passing A/C/G/T observation. '
                 'Inspect base counts and anchored indels in the site table to assess an individual locus. Reads used during polishing also contribute here.</p>')
    body += '<p>'+html.escape(meta["annotation"])+'.</p>'
    if annotations:
        body += table([{k:annotations.get(k,0) for k in ("CDS","tRNA","rRNA")}])
    body += '</div></section>'
    files = {"Stage-to-final differences":"stage_to_final.tsv", "Projected edit positions":"projected_edits.tsv", "All windows":"windows.tsv",
             "Ranked hotspots":"hotspots.tsv", "Read observations by site":"site_support.tsv", "Measurement provenance":"validation.json"}
    links = "".join(f'<li><a href="{quote(os.path.relpath(tables/path,dest))}">{label}</a></li>' for label,path in files.items())
    dest.mkdir(parents=True,exist_ok=True)
    derived = (
        ("window_complexity.tsv",complexity,("contig","start","end","full_window","edit_bases","valid_triplets","entropy_bits")),
        ("dispersion.tsv",dispersion,("comparison","full_windows","mean_edit_bases","variance_to_mean")),
        ("hotspot_motifs.tsv",motifs,("rank","contig","start","motif","window_count","window_valid_triplets","background_count","background_valid_triplets","frequency_ratio")),
        ("hotspot_context.tsv",context,("rank","contig","start","end","edit_bases","rRNA_overlap_bp","annotation_context","observed_edit_sites","mean_informative_depth","minimum_informative_depth")),
    )
    for filename, rows, fields in derived:
        out = io.StringIO()
        writer = csv.DictWriter(out,fieldnames=fields,delimiter="\t",lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        atomic_write(dest/filename,out.getvalue())
        links += f'<li><a href="{filename}">{filename}</a></li>'
    document = evidence_page("Polishing validation",iso.name,body+'<section id="validation-evidence"><h2>Evidence files</h2><ul class="evidence-links">'+links+'</ul></section>',
                             [("polish-overview","Overview"),("edit-density-atlas","Edit distribution"),("polish-context","Spatial context"),("polish-outcomes","Stage effects"),("polish-support","Read and RNA evidence"),("validation-evidence","Evidence files")], report_type="validation")
    atomic_write(dest/"summary.html",body)
    atomic_write(dest/"index.html",document)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--isolate-dir",type=Path,required=True)
    parser.add_argument("--threads",type=int,default=4)
    parser.add_argument("--out",type=Path)
    parser.add_argument("--render-only",action="store_true",help="Render retained validation tables without running analysis")
    parser.add_argument("--tables",type=Path,help="Location of retained validation tables for rendering")
    parser.add_argument("--refresh-report",action="store_true")
    args=parser.parse_args()
    iso=args.isolate_dir.resolve()
    dest=args.out or iso/"reports/polishing_validation"
    if args.threads<1:
        parser.error("threads must be positive")
    dest = dest.resolve()
    if args.tables and not args.render_only:
        parser.error("--tables requires --render-only")
    dest.parent.mkdir(parents=True,exist_ok=True)
    with (dest.parent/".polishing-validation.lock").open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if args.render_only:
            tables = (args.tables or iso/"reports/polishing_validation").resolve()
            render(iso,tables,dest)
        elif not (dest/"index.html").is_file():
            work=Path(tempfile.mkdtemp(dir=dest.parent,prefix="validation-attempt-"))
            run(iso,work,args.threads)
            work.rename(dest)
        else:
            render(iso,dest,dest)
        print(dest/"index.html")

        if args.refresh_report:
            page = iso/"reports/hybas_research/hybas_report.html"
            if page.is_file():
                original = page.read_text()
                target = quote(os.path.relpath(dest/"index.html",page.parent))
                updated = refresh_section(original,(dest/"summary.html").read_text(),target)
                atomic_write(page,updated)
                print(page)

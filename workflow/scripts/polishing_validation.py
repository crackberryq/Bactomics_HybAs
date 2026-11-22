"""Backfill final-coordinate evidence from retained stages without re-running polishing."""
import argparse
import base64
import io
import re
import csv
import fcntl
import gzip
import hashlib
import html
import json
import os
import shutil
import subprocess
import tempfile
from bisect import bisect_left
from collections import Counter, defaultdict
from pathlib import Path

from polishing_trajectory import fasta_lengths, write_table
from hybas_research_report import STYLE

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
            for contig, pos in sorted(sites):
                counts, indels = Counter(), Counter()
                for col in align.pileup(contig,pos-1,pos,truncate=True,stepper="samtools",min_mapping_quality=20,min_base_quality=20,max_depth=1000000):
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
                vals = [counts[b] for b in "ACGT"]
                rows.append(dict(zip(fields,[contig,pos,fa.fetch(contig,pos-1,pos).upper(),sum(vals),*vals,counts["del"],json.dumps(indels,sort_keys=True)])))
        status = "Final-reference Illumina pileup, mapping and base quality ≥20"
    write_table(root/"site_support.tsv",rows,fields)
    return status


def figures(stages, windows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def embed(fig, label):
        buf = io.BytesIO()
        fig.savefig(buf,format="png",dpi=125,bbox_inches="tight")
        plt.close(fig)
        image = base64.b64encode(buf.getvalue()).decode("ascii")
        return f'<img style="width:100%;height:auto" alt="{label}" src="data:image/png;base64,{image}">'

    names = [s["stage"] for s in stages]
    x = list(range(1,len(names)+1))
    short = [float(s["short_pct"]) for s in stages]
    mean = [float(s["mean_aa"]) for s in stages]
    complete = [float(s["complete_pct"]) for s in stages]
    fig = plt.figure(figsize=(10,7))
    fig.subplots_adjust(wspace=.34,hspace=.56)
    ax = fig.add_subplot(221,projection="3d")
    ax.plot(x,short,mean,color="#136d68")
    ax.scatter(x,short,mean,color="#136d68")
    ax.set(xlabel="Stage",ylabel="Short (%)",zlabel="Mean (aa)",title="Oblique trajectory")
    for panel,xx,yy,title,xlabel,ylabel in [
        (222,x,short,"Top view","Stage","Proteins under 100 aa (%)"),
        (223,x,mean,"Front view","Stage","Mean protein length (aa)"),
        (224,short,mean,"Side view","Proteins under 100 aa (%)","Mean protein length (aa)")]:
        ax = fig.add_subplot(panel)
        ax.plot(xx,yy,"o-",color="#136d68")
        ax.set(xlabel=xlabel,ylabel=ylabel,title=title)
        if panel!=224:
            ax.set_xticks(x)
        ax.grid(alpha=.2)
    view = embed(fig,"Oblique, top, front and side views of the polishing trajectory")
    view += "<p>" + " → ".join(f"{i}. {html.escape(name)}" for i,name in zip(x,names)) + "</p>"
    fig, ax = plt.subplots(figsize=(9,2.8),layout="constrained")
    ax.plot(x,complete,"o-",color="#136d68")
    ax.set_xticks(x,names)
    ax.set_ylabel("Complete BUSCO (%)")
    ax.grid(axis="y",alpha=.2)
    view += embed(fig,"BUSCO completeness by stage")
    contigs = {}
    for w in windows:
        contigs.setdefault(w["contig"],[]).append(w)
    fig,ax = plt.subplots(figsize=(9,3),layout="constrained")
    offset = 0
    for contig,rows in contigs.items():
        ax.plot([(offset+w["start"]-1)/1e6 for w in rows],[w["edit_bases"] for w in rows],lw=1)
        offset += rows[-1]["end"]
    ax.set(xlabel="Position along concatenated final contigs (Mb)",ylabel="Edited bases per window")
    return view,embed(fig,"Projected edit burden by 10 kb final-reference window")


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
    def table(rows):
        if not rows:
            return "<p>No observations.</p>"
        return "<table><tr>"+ "".join("<th>"+html.escape(k.replace("_"," "))+"</th>" for k in rows[0])+"</tr>"+ "".join("<tr>"+ "".join("<td>"+html.escape(str(v))+"</td>" for v in r.values())+"</tr>" for r in rows)+"</table>"
    views,density = figures(stages,windows)
    body = ("<h2>Polishing validation</h2><div class='panel'><h3>Each retained stage compared with the final assembly</h3>"+views+table(direct)+
            "<p>Differences cover one-to-one aligned sequence. Stage-only and final-only columns count bases. Individual alignment reports record unaligned sequence.</p>"+
            f"<h3>Spatial edit density</h3><p>{meta['projected']:,} adjacent-stage edit bases projected to the final assembly. {meta['unresolved']:,} unresolved bases were excluded. Windows span 10 kb, with a shorter terminal window per contig.</p>"+density+
            "<h3>Highest-burden windows</h3>"+table(hotspots)+
            "<p>Overlap counts rRNA bases inside each window. It does not count individual edits in rRNA.</p><p>"+html.escape(gff_status)+". "+
            html.escape(", ".join(f"{k}: {v}" for k,v in sorted(counts.items()) if k in ("rRNA","tRNA","CDS")))+
            ".</p><h3>Site-level read evidence</h3><p>"+html.escape(support)+". The site table records A/C/G/T observations, deletion spans and anchored indels. Sites falling outside unique alignments have no final-coordinate support estimate.</p></div>")
    (root/"summary.html").write_text(body)
    files = ["stage_to_final.tsv","projected_edits.tsv","windows.tsv","hotspots.tsv","site_support.tsv","validation.json"]
    links = " · ".join(f'<a href="{name}">{name}</a>' for name in files)
    (root/"index.html").write_text('<meta charset="utf-8"><title>HybAs polishing validation</title><style>'+STYLE+'</style><main>'+body+'<p>'+links+'</p></main>')


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--isolate-dir",type=Path,required=True)
    parser.add_argument("--threads",type=int,default=4)
    parser.add_argument("--out",type=Path)
    parser.add_argument("--refresh-report",action="store_true")
    args=parser.parse_args()
    iso=args.isolate_dir.resolve()
    dest=args.out or iso/"reports/polishing_validation"
    if args.threads<1:
        parser.error("threads must be positive")
    with (dest.parent/".polishing-validation.lock").open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if not (dest/"index.html").is_file():
            work=Path(tempfile.mkdtemp(dir=dest.parent,prefix="validation-attempt-"))
            run(iso,work,args.threads)
            work.rename(dest)
        print(dest/"index.html")

        if args.refresh_report:
            page = iso/"reports/hybas_research/hybas_report.html"
            if page.is_file():
                original = page.read_text()
                marker = "<!-- HybAs final-coordinate validation -->"
                if marker not in original:
                    if "</main>" not in original:
                        raise ValueError(f"Unexpected research report layout: {page}")
                    section = marker + '<section id="polishing-validation">' + (dest/"summary.html").read_text() + '<a href="../polishing_validation/index.html">Validation evidence tables</a></section>'
                    with tempfile.NamedTemporaryFile(mode="w",dir=page.parent,delete=False) as out:
                        out.write(original.replace("</main>",section+"</main>",1))
                    os.replace(out.name,page)
                print(page)

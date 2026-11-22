"""Offline vector plots for final-coordinate validation reports."""
import base64
import html
import math
import re
import statistics
from collections import Counter, defaultdict

colors={'substitution':'#3678ab','insertion':'#e08035','deletion':'#865ba4'}
def density_chart(contig, rows, pairs, events):
    rows=sorted(rows)
    end=rows[-1][1]
    left,width=195,790
    def x(pos): return left+width*(pos-1)/max(1,end-1)
    height=285+48*len(pairs)
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1030 {height}" font-family="Arial,sans-serif" fill="#21434b"><rect width="100%" height="100%" fill="white"/>']
    def label(xx,yy,txt,size=14,anchor='start'):
        parts.append(f'<text x="{xx}" y="{yy}" text-anchor="{anchor}" font-size="{size}">{html.escape(str(txt))}</text>')
    peak=max(v[2] for v in rows); scale=max(1,peak)
    label(12,32,'Density',18);label(12,55,'edited bases / window',12)
    points=[]
    for lo,hi,n,_ in rows:
        points.extend([(x(lo),142-110*n/scale), (x(hi),142-110*n/scale)])
    line=' '.join(f'{xx:.2f},{yy:.2f}' for xx,yy in points)
    parts.append(f'<polygon points="{points[0][0]},142 {line} {points[-1][0]},142" fill="#cde7e3"/><polyline points="{line}" fill="none" stroke="#087b73" stroke-width="2"/>')
    for n in sorted({0,scale//2,scale}):
        yy=142-110*n/scale;label(left-8,yy+4,n,12,'end')
        parts.append(f'<path d="M{left} {yy}h{width}" stroke="#668988" stroke-opacity=".22"/>')
    label(12,177,'rRNA overlap',15)
    known=any(r[3]!='' for r in rows)
    if known:
        parts.append(f'<rect x="{left}" y="162" width="{width}" height="18" fill="#eff4f4"/>')
        for lo,hi,_,rna in rows:
            if rna and int(rna):
                parts.append(f'<rect x="{x(lo):.2f}" y="162" width="{max(.8,x(hi)-x(lo)):.2f}" height="18" fill="#28818d"/>')
    else:label(left,177,'Annotation coordinates unavailable',13)
    for i,pair in enumerate(pairs):
        yy=213+48*i
        label(12,yy+3,pair.replace('-to-',' → '),14)
        parts.append(f'<rect x="{left}" y="{yy-16}" width="{width}" height="32" rx="4" fill="#f3f6f6"/>')
        # Display overlapping marks once at each pixel and retain exact counts in the window track.
        marks = defaultdict(set)
        for pos, kind in events[contig,pair]:
            marks[kind].add(round(x(pos),1))
        for k,(kind,color) in enumerate(colors.items()):
            segments = "".join(f"M{px} {yy+9*(k-1)-3}v6" for px in sorted(marks[kind]))
            parts.append(f'<path d="{segments}" stroke="{color}"/>')
    yy=217+48*len(pairs)
    for i in range(6):label(left+width*i/5,yy,f'{(1+(end-1)*i/5)/1e6:.3f}',12,'middle')
    label(left+width/2,yy+23,'Final assembly position (Mb)',14,'middle')
    for i,(k,c) in enumerate(colors.items()):
        xx=left+i*185
        parts.append(f'<rect x="{xx}" y="{yy+40}" width="12" height="6" fill="{c}"/>')
        label(xx+18,yy+47,k.capitalize(),13)
    parts.append('</svg>')
    return base64.b64encode(''.join(parts).encode()).decode()


def density_atlas(windows, edits, stages):
    grouped = defaultdict(list)
    for row in windows:
        grouped[row["contig"]].append((int(row["start"]), int(row["end"]), int(row["edit_bases"]), row["rRNA_overlap_bp"]))
    names = [r["stage"] for r in stages]
    pairs = [a+"-to-"+b for a,b in zip(names,names[1:])]
    events = defaultdict(list)
    for row in edits:
        if row["status"]=="mapped":
            events[row["final_contig"],row["transition"]].append((int(row["final_pos"]),row["event"]))
    charts = []
    empty = 0
    for contig,rows in sorted(grouped.items(),key=lambda x:-sum(r[2] for r in x[1])):
        lo,hi,peak,_ = max(rows,key=lambda r:r[2])
        if not peak:
            empty += 1
            continue
        caption = f'{html.escape(contig)} · highest window {lo:,}–{hi:,} bp ({peak:,} edited bases)'
        fig = density_chart(contig,rows,pairs,events)
        img = f'<img style="display:block;width:100%;height:auto" alt="Density, rRNA-positive windows and edits by stage" src="data:image/svg+xml;base64,{fig}">'
        if charts:
            charts.append(f'<details><summary>{caption}</summary>{img}</details>')
        else:
            charts.append(f'<h4>{caption}</h4>{img}')
    if empty:
        charts.append(f'<p>{empty:,} contig(s) had zero projected edits.</p>')
    return "".join(charts)


def projections(stages):
    names = [r["stage"] for r in stages]
    order = list(range(1,len(names)+1))
    short = [float(r["short_pct"]) for r in stages]
    mean = [float(r["mean_aa"]) for r in stages]
    busco = [float(r["complete_pct"]) for r in stages]
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 930 486" font-family="system-ui,sans-serif" fill="#183944"><rect width="100%" height="100%" fill="white"/>']
    def label(x,y,value,size=13,anchor="start"):
        parts.append(f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" font-size="{size}">{html.escape(str(value))}</text>')
    for i,(view,xx,yy,xlabel,ylabel) in enumerate([
        ("Top view",order,short,"Stage","Short proteins (%)"),
        ("Front view",order,mean,"Stage","Mean protein length (aa)"),
        ("Side view",short,mean,"Proteins under 100 aa (%)","Mean protein length (aa)"),
        ("Complete BUSCO",order,busco,"Stage","Complete markers (%)"),
    ]):
        left,top = 68+(i%2)*460,60+(i//2)*240
        label(left,top-35,view,19)
        label(left,top-14,ylabel,13)
        def scale(v,values,span):
            delta=max(values)-min(values)
            return ((v-min(values))/delta if delta else .5)*span
        for tick in sorted({min(yy),(min(yy)+max(yy))/2,max(yy)}):
            y=top+120-scale(tick,yy,115)
            label(left-8,y+4,f"{tick:.1f}",12,"end")
            parts.append(f'<path d="M{left} {y:.1f}h337" stroke="#d8e4e2"/>')
        for tick in order if xlabel=="Stage" else sorted({min(xx),(min(xx)+max(xx))/2,max(xx)}):
            label(left+scale(tick,xx,335),top+146,str(tick) if xlabel=="Stage" else f"{tick:.1f}",12,"middle")
        label(left+167,top+169,xlabel,13,"middle")
        coords=[(left+scale(x,xx,335),top+120-scale(y,yy,115)) for x,y in zip(xx,yy)]
        path=" ".join(f"{x:.1f},{y:.1f}" for x,y in coords)
        parts.append(f'<polyline points="{path}" fill="none" stroke="#117d73" stroke-width="2.5"/>')
        for j,(x,y) in enumerate(coords):
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="#117d73"/>')
            if view=="Side view":
                label(x+7,y+[0,12,-4,8,14][j%5],j+1,11)
    parts.append('</svg>')
    data=base64.b64encode("".join(parts).encode()).decode()
    return (f'<img style="width:100%" alt="Top, front and side projections of the trajectory plus BUSCO by stage" src="data:image/svg+xml;base64,{data}">'
            '<p>'+" → ".join(f'{n} {html.escape(stage)}' for n,stage in enumerate(names,1))+'</p>')


def spatial_context(reference, windows, hotspots, edits, support, pairs, size):
    """Use identical final-coordinate bins for retained and newly measured evidence."""
    counts_by_bin = defaultdict(Counter)
    for row in edits:
        if row["status"] == "mapped":
            pos = int(row["final_pos"])
            contig = row["final_contig"]
            pair = row["transition"]
            if pair not in pairs or not 1 <= pos <= len(reference[contig]):
                raise ValueError("Projected edit outside recorded stages or final sequence")
            counts_by_bin[contig, 1 + ((pos - 1) // size) * size][pair] += 1
    ranked = {(r["contig"], int(r["start"])): i for i, r in enumerate(hotspots, 1)}
    required = {(c, i+1) for c,seq in reference.items() for i in range(0, len(seq), size)}
    seen = set()
    profiles, complexity, pooled, by_pair = {}, [], Counter(), defaultdict(list)
    triplets = re.compile(r"(?=([ACGT]{3}))")
    for row in windows:
        contig, start, end = row["contig"], int(row["start"]), int(row["end"])
        key = contig, start
        burden = int(row["edit_bases"])
        if (key not in required or key in seen or end != min(start+size-1, len(reference[contig]))
                or sum(counts_by_bin[key].values()) != burden):
            raise ValueError("Window boundaries or totals differ from projected edit table")
        seen.add(key)
        motif_counts = Counter(triplets.findall(reference[contig][start-1:end]))
        total = sum(motif_counts.values())
        entropy = -sum((n/total)*math.log2(n/total) for n in motif_counts.values())
        complete = end-start+1 == size
        complexity.append(dict(contig=contig, start=start, end=end, full_window=int(complete), edit_bases=burden,
                               valid_triplets=total, entropy_bits=round(max(0, entropy), 6) if total else ""))
        if key in ranked:
            profiles[key] = motif_counts
        else:
            pooled.update(motif_counts)
        if complete:
            for pair in pairs:
                by_pair[pair].append(counts_by_bin[key][pair])
            by_pair["All transitions"].append(burden)
    if seen != required or len(profiles) != len(ranked):
        raise ValueError("Spatial evidence does not cover the retained final sequence")
    dispersion = []
    for name in [*pairs, "All transitions"]:
        bins = by_pair[name]
        mean = statistics.mean(bins) if bins else None
        vmr = statistics.variance(bins)/mean if len(bins)>1 and mean else None
        dispersion.append(dict(comparison=name, full_windows=len(bins),
                               mean_edit_bases=round(mean,6) if mean is not None else "",
                               variance_to_mean=round(vmr,6) if vmr is not None else ""))
    edited_sites = {(r["final_contig"],int(r["final_pos"])) for r in edits if r["status"] == "mapped"}
    local_depth = defaultdict(list)
    for row in support:
        pos = int(row["position"])
        if (row["contig"],pos) in edited_sites:
            local_depth[row["contig"],1+((pos-1)//size)*size].append(int(row["informative_bases"]))
    motifs, context = [], []
    bg_total = sum(pooled.values())
    for row in hotspots:
        key = row["contig"],int(row["start"])
        motif_counts = profiles[key]
        total = sum(motif_counts.values())
        rrna = row["rRNA_overlap_bp"]
        label = "Annotation unavailable" if rrna == "" else ("Overlaps rRNA" if int(rrna) else "Outside annotated rRNA")
        depths = local_depth[key]
        context.append(dict(rank=ranked[key], **row, annotation_context=label, observed_edit_sites=len(depths),
                            mean_informative_depth=round(statistics.mean(depths),2) if depths else "",
                            minimum_informative_depth=min(depths) if depths else ""))
        if bg_total:
            for motif, count in sorted(motif_counts.items()):
                fold = ((count+.5)/(total+32))/((pooled[motif]+.5)/(bg_total+32))
                if fold > 1:
                    motifs.append(dict(rank=ranked[key], contig=key[0], start=key[1], motif=motif, window_count=count,
                                       window_valid_triplets=total, background_count=pooled[motif], background_valid_triplets=bg_total,
                                       frequency_ratio=round(fold,6)))
    return complexity, dispersion, motifs, context


def complexity_plot(rows):
    values=[(float(r["entropy_bits"]),int(r["edit_bases"])) for r in rows if r["full_window"] and r["entropy_bits"] != ""]
    if not values:
        return "<p>No full windows contained measurable triplets.</p>"
    low, high = min(x for x,_ in values), max(x for x,_ in values)
    peak = max(y for _,y in values)
    def sx(v): return 86+722*((v-low)/(high-low) if high>low else .5)
    def sy(v): return 232-190*v/max(1,peak)
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 864 300" font-family="system-ui,sans-serif" fill="#183944"><rect width="100%" height="100%" fill="white"/><text x="16" y="20">Projected edited bases / window</text>']
    for t in sorted({0,peak//2,peak}):
        y = sy(t)
        parts.append(f'<path d="M84 {y:.1f}h728" stroke="#dae6e4"/><text x="72" y="{y+5:.1f}" text-anchor="end">{t}</text>')
    for t in sorted({low,(low+high)/2,high}):
        parts.append(f'<text x="{sx(t):.2f}" y="256" text-anchor="middle">{t:.3f}</text>')
    for x,y in values:
        parts.append(f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="3.4" fill="#117c78" fill-opacity=".45"/>')
    parts.append('<text x="448" y="285" text-anchor="middle">3-mer Shannon entropy (bits)</text></svg>')
    data=base64.b64encode("".join(parts).encode()).decode()
    return f'<img style="width:100%" alt="Genome-wide full windows by sequence complexity and projected edit burden" src="data:image/svg+xml;base64,{data}">'

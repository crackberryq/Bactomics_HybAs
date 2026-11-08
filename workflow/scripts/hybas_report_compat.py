"""Compatibility readers for historical HybAs evidence; never execute tools."""
from __future__ import annotations
import csv
import re
from pathlib import Path


def summary_tsv(path):
    if not path.is_file():
        return {}
    with path.open(encoding="utf8", errors="replace") as f:
        return {row["metric"]: row["value"] for row in csv.DictReader(f, delimiter="\t")
                if row.get("metric") and row.get("value") is not None}


def quast(path):
    """Parse QUAST metric-per-row report.txt, not its transposed representation."""
    if not path.is_file():
        return {}
    fields = {"# contigs": "assembly_contig_count_fasta", "Total length": "assembly_bp",
              "N50": "quast_n50", "L50": "quast_l50", "GC (%)": "quast_gc_percent",
              "Largest contig": "quast_largest_contig"}
    result = {}
    with path.open(encoding="utf8", errors="replace") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                parts = re.split(r"\s{2,}", line.strip())
            if len(parts) < 2:
                continue
            key = parts[0].strip()
            if key in fields:
                result[fields[key]] = parts[1].strip()
    return result


def busco(path):
    if not path.is_file():
        return {}
    text = path.read_text(encoding="utf8", errors="replace")
    match = re.search(r"C:\s*([\d.]+)%.*?S:\s*([\d.]+)%.*?D:\s*([\d.]+)%.*?F:\s*([\d.]+)%.*?M:\s*([\d.]+)%", text, re.S)
    if not match:
        return {}
    return dict(zip(("busco_complete_percent", "busco_single_copy_percent",
                     "busco_duplicated_percent", "busco_fragmented_percent",
                     "busco_missing_percent"), match.groups()))


def busco_path(isolate):
    preferred = isolate / "reports/busco"
    matches = list(preferred.glob("**/short_summary*.txt")) if preferred.is_dir() else []
    # Only treat a final BUSCO result as final. Do not silently use intermediate results.
    matches = [p for p in matches if "final" in str(p).lower() or "busco_" + isolate.name in str(p)]
    if matches:
        return sorted(matches, key=lambda p: (len(p.parts), str(p)))[0]
    for root in (isolate / "reports/cds_audit/busco/final",):
        if root.is_dir():
            options = sorted(root.glob("**/short_summary*.txt"))
            if options:
                return options[0]
    return None


def stage_quast(isolate):
    """Use one preferred stage-QC collection; avoid counting duplicate audits."""
    roots = [isolate / "reports/FINAL_AUDIT_RESULTS/QC_ASSEMBLY/qc_steps",
             isolate / "reports/intermediate_qc_FULL",
             isolate / "reports/full_audit/qc_steps"]
    for root in roots:
        if not root.is_dir():
            continue
        stages = []
        for d in sorted(root.iterdir()):
            report = d / "report.txt"
            if d.is_dir() and report.is_file():
                m = quast(report)
                if m:
                    stages.append({"stage": d.name,
                                   "contigs": m.get("assembly_contig_count_fasta", "—"),
                                   "assembly_bp": m.get("assembly_bp", "—"),
                                   "n50": m.get("quast_n50", "—"),
                                   "gc": m.get("quast_gc_percent", "—")})
        if stages:
            return stages, root
    return [], None


def evidence(isolate, summary):
    sources = {}
    def add(key, data, path):
        if data:
            sources[key] = str(path.relative_to(isolate))
            for k, v in data.items():
                if summary.get(k) in (None, "", "NA", "Not recorded"):
                    summary[k] = v
    path = isolate / "metadata/run_summary.tsv"
    add("run_summary_tsv", summary_tsv(path), path)
    path = isolate / "reports/quast/report.txt"
    add("quast", quast(path), path)
    bp = busco_path(isolate)
    if bp:
        add("busco", busco(bp), bp)
    stages, root = stage_quast(isolate)
    if root:
        sources["stage_quast"] = str(root.relative_to(isolate))
    return summary, stages, sources

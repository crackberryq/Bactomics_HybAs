#!/usr/bin/env python3
"""Standalone, offline HybAs Research Report. Never schedules analysis tools."""
from __future__ import annotations

import argparse
import csv
import html
import json
from datetime import datetime, timezone
from pathlib import Path

from aware_read_audit import audit, audit_ont, decision


STYLE = """
:root{--navy:#112d35;--ink:#1b2d40;--teal:#087c78;--paper:#f4f6f7;--line:#dce3e8;--muted:#5e6c79}
*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:14px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
.app{display:grid;grid-template-columns:230px minmax(0,1fr);min-height:100vh}
aside{background:var(--navy);color:#e1e9ec;padding:28px 18px;position:sticky;top:0;height:100vh}
.brand{font-size:27px;font-weight:750;color:white;margin-bottom:6px}.sub{font-size:10px;letter-spacing:1.5px;color:#9fc6c5;margin-bottom:32px}
aside a{display:block;color:#c7d8dc;text-decoration:none;padding:12px;border-radius:7px}aside a:hover{background:#25444d;color:white}
main{padding:30px 38px;max-width:1500px;width:100%}
.hero{background:linear-gradient(115deg,#142d43,#173e52);border-radius:17px;color:white;padding:35px;margin-bottom:24px}
.hero h1{font-size:34px;margin:0 0 8px}.hero p{margin:0;color:#c9e3e6}
.grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}
.stat,.panel{background:white;border:1px solid var(--line);border-radius:12px;padding:20px;margin-bottom:18px}
.stat small{color:var(--muted);display:block}.stat strong{font-size:25px;display:block;overflow-wrap:anywhere}
section{margin:34px 0}h2{font-size:23px}h3{font-size:16px}p.note{color:var(--muted)}
table{width:100%;border-collapse:collapse;font-size:12px}th{text-align:left;background:#edf3f2;color:#38565c}
td,th{border-bottom:1px solid var(--line);padding:10px;vertical-align:top;overflow-wrap:anywhere}
.table-wrap{overflow:auto;max-height:650px}tbody tr:hover{background:#f5f9f8}
code{overflow-wrap:anywhere}.badge{display:inline-block;background:#e6f2ef;color:#087c78;padding:4px 10px;border-radius:5px;font-size:12px}
.warn{background:#fff5df;color:#80540c;padding:12px;border-radius:8px}
input{width:100%;padding:10px;border:1px solid var(--line);border-radius:7px;margin-bottom:10px}
a{color:var(--teal)}footer{font-size:11px;color:var(--muted);padding:30px 0}
@media(max-width:850px){.app{display:block}aside{position:static;height:auto}main{padding:16px}.grid{grid-template-columns:repeat(2,1fr)}}
"""


def read_json(path):
    try:
        return json.loads(path.read_text(encoding="utf8"))
    except (OSError, ValueError):
        return {}


def esc(x):
    return html.escape(str(x if x is not None else "Not recorded"), quote=True)


def table(columns, rows):
    head = "".join("<th>" + esc(k) + "</th>" for k in columns)
    body = "".join("<tr>" + "".join("<td>" + esc(row.get(k, "")) + "</td>" for k in columns) + "</tr>"
                   for row in rows)
    return "<div class='table-wrap'><table><thead><tr>" + head + "</tr></thead><tbody>" + body + "</tbody></table></div>"


def key_values(data, keys):
    return table(["Field", "Value"], [{"Field": key.replace("_", " "), "Value": data.get(key, "Not recorded")}
                                     for key in keys])


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def generate(isolate, output_dir, with_audit):
    output_dir.mkdir(parents=True, exist_ok=True)
    dec = decision(isolate / "metadata/kraken_decision.tsv")
    summary = read_json(isolate / "metadata/run_summary.json")
    execution = read_json(isolate / "metadata/polish_execution.json")
    result = audit(isolate, output_dir, dec.get("selected_target_taxid", "")) if with_audit else {
        "status": "not_requested", "note": "Use --audit-reads to reconcile merged and cleaned FASTQs."}
    ont = audit_ont(isolate, output_dir, dec.get("selected_target_taxid", "")) if with_audit else {"status": "not_requested"}
    (output_dir / "aware_ont_filter_audit.json").write_text(json.dumps(ont, indent=2) + "\n", encoding="utf8")
    if ont["status"] == "completed":
        write_tsv(output_dir / "aware_ont_removed_taxa.tsv", ont["by_taxon"],
                  ["taxid", "name", "rank", "category", "input_reads", "retained_reads", "removed_reads"])
        write_tsv(output_dir / "aware_ont_removed_categories.tsv", ont["by_category"],
                  ["category", "input_reads", "retained_reads", "removed_reads"])
    result["declared_illumina_action"] = dec.get("illumina_kraken_decision", "Not recorded")
    result["declared_ont_action"] = dec.get("ont_kraken_decision", "Not recorded")
    if result["status"] == "completed":
        t = result["totals"]
        observed_passthrough = result["integrity_ok"] and t["removed_pairs"] == 0
        result["action_consistent"] = (
            observed_passthrough if result["declared_illumina_action"] == "passthrough"
            else result["integrity_ok"] if result["declared_illumina_action"] == "whitelist_filter"
            else False)
        result["action_verification_scope"] = (
            "Observed read retention only; descendant whitelist membership is not independently verified.")
        write_tsv(output_dir / "aware_removed_taxa.tsv", result["by_taxon"],
                  ["taxid", "name", "rank", "category", "input_pairs", "retained_pairs",
                   "removed_pairs", "partial_pairs"])
        write_tsv(output_dir / "aware_removed_categories.tsv", result["by_category"],
                  ["category", "input_pairs", "retained_pairs", "removed_pairs", "partial_pairs"])
    (output_dir / "aware_filter_audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    overview = {
        "Run mode": summary.get("run_mode", "Not recorded"),
        "Assembly size": summary.get("assembly_bp", "Not recorded"),
        "N50": summary.get("quast_n50", "Not recorded"),
        "Polishing": ", ".join(execution.get("executed_stages", [])) or "None recorded"
    }
    stats = "".join("<div class='stat'><small>" + esc(k) + "</small><strong>" + esc(v) + "</strong></div>"
                    for k, v in overview.items())
    keys = ["kraken_decision_mode", "selected_target_taxid", "selected_target_name",
            "selected_target_rank", "selected_target_source", "target_selection_reason",
            "illumina_selected_target_percent", "ont_selected_target_percent",
            "illumina_kraken_decision", "illumina_kraken_decision_reason",
            "ont_kraken_decision", "ont_kraken_decision_reason",
            "selected_target_percent_for_decision", "decision_confidence", "awareness_state",
            "runner_up_taxid", "runner_up_name", "runner_up_percent",
            "selected_vs_runner_up_delta_percent", "platform_agreement_score",
            "platform_agreement_reason", "kraken_include_children"]
    audit_html = "<p class='warn'>Read audit: " + esc(result["status"]) + "</p>"
    if result["status"] == "completed":
        totals = result["totals"]
        audit_html = ("<p class='badge'>Paired-read reconciliation: " +
                      ("PASS" if result["integrity_ok"] else "ATTENTION") + "</p>" +
                      key_values(totals, list(totals)) +
                      "<h3>Removed reads by taxonomic category</h3>" +
                      table(["category", "input_pairs", "retained_pairs", "removed_pairs"],
                            result["by_category"]) +
                      "<h3>Taxon-level detail</h3><input id='taxon-search' placeholder='Search taxon, TaxID, category...' aria-label='Filter taxa'>" +
                      "<div id='taxon-table'>" +
                      table(["taxid", "name", "rank", "category", "input_pairs", "retained_pairs", "removed_pairs"],
                            result["by_taxon"][:500]) + "</div>" +
                      "<p class='note'>First 500 taxa displayed. Complete data are exported to aware_removed_taxa.tsv.</p>" +
                      "<p class='note'>Action consistency: " + esc(result["action_consistent"]) +
                      ". This is not an independent verification of whitelist lineage membership.</p>")
    elif result.get("missing_files"):
        audit_html += "<p class='note'>Missing inputs: " + esc(", ".join(result["missing_files"])) + "</p>"
    ont_html = "<p class='warn'>ONT audit: " + esc(ont["status"]) + "</p>"
    if ont["status"] == "completed":
        ont_html = ("<p class='badge'>ONT read reconciliation: " +
                    ("PASS" if ont["integrity_ok"] else "ATTENTION") + "</p>" +
                    key_values(ont["totals"], list(ont["totals"])) +
                    "<h3>ONT removed reads by category</h3>" +
                    table(["category", "input_reads", "retained_reads", "removed_reads"], ont["by_category"]) +
                    "<h3>ONT taxonomic assignments</h3>" +
                    table(["taxid", "name", "category", "input_reads", "retained_reads", "removed_reads"],
                          ont["by_taxon"][:500]))
    stages = execution.get("stages", [])
    stage_rows = [{"Stage": s.get("stage"), "FASTA bases": s.get("bases"),
                   "Contigs": s.get("contigs"), "SHA256": s.get("sha256")} for s in stages]
    generated = datetime.now(timezone.utc).isoformat()
    body = (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        "<title>HybAs · " + esc(isolate.name) + "</title><style>" + STYLE + "</style></head><body>"
        "<div class='app'><aside><div class='brand'>HybAs</div><div class='sub'>BACTOMICS · RESEARCH EVIDENCE</div>"
        "<a href='#overview'>01 · Overview</a><a href='#aware'>02 · AWARE decisions</a>"
        "<a href='#removed'>03 · Removed reads</a><a href='#assembly'>04 · Assembly & polishing</a>"
        "<a href='#provenance'>05 · Evidence files</a></aside><main>"
        "<div class='hero'><h1>" + esc(isolate.name) + "</h1><p>HybAs Research Report · local, read-only evidence dashboard</p></div>"
        "<section id='overview'><h2>Overview</h2><div class='grid'>" + stats + "</div></section>"
        "<section id='aware'><h2>AWARE decision explorer</h2><div class='panel'>"
        "<p class='note'>Target election, per-platform read actions and conservative confidence are separate outputs. "
        "A low confidence label does not override a platform's passthrough.</p>" +
        key_values(dec, keys) + "</div></section>"
        "<section id='removed'><h2>Removed-read origins</h2><div class='panel'>" + audit_html +
        "<h3>ONT filtering</h3>" + ont_html +
        "<p class='note'>Kraken classifications are taxonomic assignments, not confirmed biological origin. "
        "Paired reads are counted once per normalized pair identifier. No read sequences or IDs are exported.</p>"
        "</div></section>"
        "<section id='assembly'><h2>Assembly and polishing</h2><div class='panel'>" +
        key_values(execution, ["requested_polish_mode", "run_mode", "planned_stages",
                               "executed_stages", "skipped_stages", "final_source_stage",
                               "final_source_path", "final_source_sha256"]) +
        "<h3>Recorded assembly stages</h3>" +
        table(["Stage", "FASTA bases", "Contigs", "SHA256"], stage_rows) + "</div></section>"
        "<section id='provenance'><h2>Evidence and exports</h2><div class='panel'>"
        "<p><a href='aware_filter_audit.json'>Filtering audit JSON</a> · "
        "<a href='aware_removed_taxa.tsv'>Taxon-level TSV</a> · "
        "<a href='aware_removed_categories.tsv'>Category TSV</a> · " +
        "<a href='aware_ont_filter_audit.json'>ONT audit JSON</a> · " +
        "<a href='aware_ont_removed_taxa.tsv'>ONT taxon TSV</a></p>"
        "<p class='note'>TSV exports exist only when a complete read audit was possible. "
        "Source files: metadata/kraken_decision.tsv, metadata/run_summary.json, "
        "metadata/polish_execution.json, reports/kraken2, work/inputs_merged, work/inputs_cleaned.</p>"
        "</div></section><footer>Generated " + esc(generated) +
        " UTC · HybAs report uses recorded evidence only · no classification or assembly was executed</footer>"
        "</main></div><script>const q=document.getElementById('taxon-search');"
        "if(q){q.addEventListener('input',()=>{document.querySelectorAll('#taxon-table tbody tr')"
        ".forEach(r=>{r.hidden=!r.textContent.toLowerCase().includes(q.value.toLowerCase())})})}"
        "</script></body></html>"
    )
    (output_dir / "hybas_report.html").write_text(body, encoding="utf8")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base-dir", required=True, type=Path)
    p.add_argument("--isolate", required=True)
    p.add_argument("--outdir", type=Path)
    p.add_argument("--audit-reads", action="store_true",
                   help="Disk-backed comparison of merged/cleaned FASTQ IDs and saved Kraken assignments")
    args = p.parse_args()
    if not args.isolate or args.isolate in (".", "..") or "/" in args.isolate or "\\" in args.isolate:
        p.error("Invalid isolate identifier")
    iso = args.base_dir / args.isolate
    if not iso.is_dir():
        p.error("Isolate directory does not exist: " + str(iso))
    output = args.outdir or iso / "reports/hybas_research"
    result = generate(iso, output, args.audit_reads)
    print("Report:", output / "hybas_report.html")
    print("Read audit:", result["status"])


if __name__ == "__main__":
    main()

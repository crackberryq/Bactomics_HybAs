#!/usr/bin/env python3
"""Standalone, offline HybAs Research Report. Never schedules analysis tools."""
from __future__ import annotations

import argparse
import csv
import html
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

from aware_read_audit import audit, audit_long, decision
from hybas_report_compat import evidence

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from report_theme import STYLE


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


def generate(isolate, output_dir, with_audit, trajectory_dir=None, validation_dir=None, audit_from=None):
    trajectory_dir = Path(trajectory_dir or isolate / "reports/polishing_trajectory")
    validation_dir = Path(validation_dir or isolate / "reports/polishing_validation")
    def href(path):
        return esc(quote(os.path.relpath(path, output_dir)))
    def saved_audit(filename):
        path = Path(audit_from) / filename
        if not path.is_file():
            return {"status": "unavailable", "note": "No retained audit record"}
        data = json.loads(path.read_text(encoding="utf8"))
        if not isinstance(data, dict) or not data.get("status"):
            raise ValueError(f"Invalid retained filtering audit: {path}")
        return data
    output_dir.mkdir(parents=True, exist_ok=True)
    dec = decision(isolate / "metadata/kraken_decision.tsv")
    summary = read_json(isolate / "metadata/run_summary.json")
    summary, historical_stages, evidence_sources = evidence(isolate, summary)
    execution = read_json(isolate / "metadata/polish_execution.json")
    platform = "pacbio" if summary.get("long_read_platform") == "pacbio" else "ont"
    platform_label = "PacBio" if platform == "pacbio" else "ONT"
    if audit_from:
        result = saved_audit("aware_filter_audit.json")
        ont = saved_audit(f"aware_{platform}_filter_audit.json")
    elif with_audit:
        target = dec.get("selected_target_taxid", "")
        result = audit(isolate, output_dir, target)
        ont = audit_long(isolate, output_dir, target, platform)
    else:
        result = {"status": "not_requested", "note": "Use --audit-reads to reconcile merged and cleaned FASTQs."}
        ont = {"status": "not_requested"}
    (output_dir / f"aware_{platform}_filter_audit.json").write_text(json.dumps(ont, indent=2) + "\n", encoding="utf8")
    if ont["status"] == "completed":
        write_tsv(output_dir / f"aware_{platform}_removed_taxa.tsv", ont["by_taxon"],
                  ["taxid", "name", "rank", "category", "input_reads", "retained_reads", "removed_reads"])
        write_tsv(output_dir / f"aware_{platform}_removed_categories.tsv", ont["by_category"],
                  ["category", "input_reads", "retained_reads", "removed_reads"])
    if not audit_from:
        result["declared_illumina_action"] = dec.get("illumina_kraken_decision", "Not recorded")
        result[f"declared_{platform}_action"] = dec.get(f"{platform}_kraken_decision", "Not recorded")
        if result["status"] == "completed":
            t = result["totals"]
            observed_passthrough = result["integrity_ok"] and t["removed_pairs"] == 0
            result["action_consistent"] = (
                observed_passthrough if result["declared_illumina_action"] == "passthrough"
                else result["integrity_ok"] if result["declared_illumina_action"] == "whitelist_filter"
                else False)
            result["action_verification_scope"] = (
                "Observed read retention only; descendant whitelist membership is not independently verified.")
    if result["status"] == "completed":
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
        "BUSCO complete": (str(summary["busco_complete_percent"]) + "%") if summary.get("busco_complete_percent") not in (None, "", "NA", "Not recorded") else "Not recorded"
    }
    stats = "".join("<div class='stat'><small>" + esc(k) + "</small><strong>" + esc(v) + "</strong></div>"
                    for k, v in overview.items())
    keys = ["kraken_decision_mode", "selected_target_taxid", "selected_target_name",
            "selected_target_rank", "selected_target_source", "target_selection_reason",
            "illumina_selected_target_percent", f"{platform}_selected_target_percent",
            "illumina_kraken_decision", "illumina_kraken_decision_reason",
            f"{platform}_kraken_decision", f"{platform}_kraken_decision_reason",
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
    ont_html = f"<p class='warn'>{platform_label} audit: " + esc(ont["status"]) + "</p>"
    if ont["status"] == "completed":
        ont_html = (f"<p class='badge'>{platform_label} read reconciliation: " +
                    ("PASS" if ont["integrity_ok"] else "ATTENTION") + "</p>" +
                    key_values(ont["totals"], list(ont["totals"])) +
                    f"<h3>{platform_label} removed reads by category</h3>" +
                    table(["category", "input_reads", "retained_reads", "removed_reads"], ont["by_category"]) +
                    f"<h3>{platform_label} taxonomic assignments</h3>" +
                    table(["taxid", "name", "category", "input_reads", "retained_reads", "removed_reads"],
                          ont["by_taxon"][:500]))
    stages = execution.get("stages", [])
    stage_rows = [{"Stage": s.get("stage"), "FASTA bases": s.get("bases"),
                   "Contigs": s.get("contigs"), "SHA256": s.get("sha256")} for s in stages]
    historical_rows = [{"Stage": x["stage"], "Contigs": x["contigs"],
                        "Assembly (bp)": x["assembly_bp"], "N50 (bp)": x["n50"],
                        "GC (%)": x["gc"]} for x in historical_stages]
    def numeric(value):
        try:
            return float(str(value).replace(",", ""))
        except (TypeError, ValueError):
            return None
    numeric_n50 = [numeric(x["n50"]) for x in historical_stages]
    valid_n50 = [v for v in numeric_n50 if v is not None and v >= 0]
    scale = max(valid_n50) if valid_n50 else 0
    stage_graph = "".join(
        "<div class='stage-bar-row'><span>" + esc(x["stage"]) + "</span>"
        "<div class='stage-track'><div class='stage-fill' style='width:" +
        str(round(100 * v / scale, 2) if scale and v is not None and v >= 0 else 0) +
        "%'></div></div><strong>" + esc(x["n50"]) + "</strong></div>"
        for x, v in zip(historical_stages, numeric_n50)
    )
    stage_section = ("<h3>Assembly-stage N50 comparison</h3><p class='note'>"
                     "Historical QUAST evidence; stages are observations, not proof that a stage executed "
                     "in this report's current environment.</p><div class='stage-graph'>" +
                     stage_graph + "</div>" +
                     table(["Stage", "Contigs", "Assembly (bp)", "N50 (bp)", "GC (%)"],
                           historical_rows)) if historical_rows else (
                     "")
    trajectory_page = trajectory_dir / "index.html"
    if trajectory_page.is_file():
        stage_section += (f'<p><a href="{href(trajectory_page)}">Explore the polishing trajectory</a> · '
                          'Compare BUSCO markers, short proteins and mean protein length across retained stages.</p>')
    validation_fragment = validation_dir / "summary.html"
    validation_section = ""
    if validation_fragment.is_file():
        validation_section = ('<!-- HybAs final-coordinate validation --><section id="polishing-validation">' +
                              validation_fragment.read_text(encoding="utf8"))
        if (validation_dir / "index.html").is_file():
            validation_section += f'<p><a href="{href(validation_dir / "index.html")}">Open polishing validation and evidence tables</a></p>'
        validation_section += '</section><!-- End HybAs final-coordinate validation -->'
    elif (validation_dir / "index.html").is_file():
        validation_section = (f'<section id="polishing-validation"><h2>Polishing validation</h2><p>'
                              f'<a href="{href(validation_dir / "index.html")}">Open polishing validation</a></p></section>')
    provenance_rows = [{"Evidence": k, "Source": v} for k, v in sorted(evidence_sources.items())]
    if not dec:
        decision_notice = ("<p class='warn'>Historical run: AWARE decision metadata was not found. "
                           "Kraken classifications may still be available, but a target decision, "
                           "confidence or whitelist action cannot be reconstructed reliably.</p>")
    else:
        decision_notice = ""
    strategy = summary.get("input_strategy", {})
    strategy_html = "<details><summary>Input libraries and assembly choice</summary><pre>" + esc(json.dumps(strategy, indent=2)) + "</pre></details>"
    model_html = "<h3>Medaka model provenance</h3>" + key_values(summary, [
        "medaka_model_requested", "medaka_model_effective", "medaka_model_sha256",
        "medaka_model_selection", "medaka_model_selection_policy", "medaka_model_selection_reason",
        "medaka_model_fallback_default", "medaka_model_fallback_medaka_version",
        "medaka_model_fallback_source_commit", "medaka_input_metadata_status"])
    if summary.get('medaka_model_selection') == 'fallback':
        model_html += ("<p class='warn'>Documented default fallback was used; this was not a match to the read metadata "
                       "and does not establish scientific suitability.</p>")
    model_html += ("<p class='note'>An explicitly selected model with missing or mixed read metadata requires researcher review. "
                   "Successful execution does not establish improved accuracy. Missing legacy model provenance is not inferred.</p>")
    compute_file = isolate.parent / "provenance/runtime/computational_report.html"
    compute_html = (f"<p><a href='{href(compute_file)}'>Computational resources and invocation history</a></p>"
                    if compute_file.is_file() else "<p class='note'>Computational resources are recorded at workflow completion.</p>")
    exports = [("aware_filter_audit.json", "Illumina audit JSON"),
               (f"aware_{platform}_filter_audit.json", f"{platform_label} audit JSON")]
    for prefix, label, audit_record in [("aware", "Illumina", result), (f"aware_{platform}", platform_label, ont)]:
        if audit_record["status"] == "completed":
            exports.extend([(f"{prefix}_removed_taxa.tsv", f"{label} taxa TSV"),
                            (f"{prefix}_removed_categories.tsv", f"{label} categories TSV")])
    file_links = "<ul class='evidence-links'>" + "".join(f"<li><a href='{name}'>{esc(label)}</a></li>" for name,label in exports) + "</ul>"
    nav_entries = [("overview","Overview"),("assembly","Assembly & polishing")]
    if validation_section:
        nav_entries.append(("polishing-validation","Polishing validation"))
    nav_entries += [("aware","AWARE decisions"),("removed","Removed reads"),("provenance","Evidence files")]
    navigation = "".join(f"<a href='#{anchor}'>{i:02d} · {esc(label)}</a>" for i,(anchor,label) in enumerate(nav_entries,1))
    generated = datetime.now(timezone.utc).isoformat()
    body = (
        "<!doctype html><html lang='en' data-hybas-report='research'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        "<title>HybAs · " + esc(isolate.name) + "</title><style>" + STYLE + "</style></head><body>"
        "<div class='app'><aside><div class='brand'>HybAs</div><div class='sub'>BACTOMICS · RESEARCH EVIDENCE</div>"
        "<nav aria-label='Report sections'>" + navigation + "</nav></aside><main>"
        "<div class='hero'><h1>" + esc(isolate.name) + "</h1><p>HybAs 10 — Experimental · AWARE+ research report · input modes and polishing require validation</p></div>"
        "<section id='overview'><h2>Overview</h2><div class='grid'>" + stats + "</div>"
        "<p class='warn'><strong>Analytical suitability: not assessed.</strong> Recorded measurements and successful processing "
        "do not establish genome correctness or suitability for downstream negative findings.</p></section>"
        "<section id='assembly'><h2>Assembly and polishing</h2><div class='panel'>" + stage_section +
        key_values(execution, ["requested_polish_mode", "executed_stages", "skipped_stages", "final_source_stage"]) +
        "<details><summary>Stage files and checksums</summary>" +
        table(["Stage", "FASTA bases", "Contigs", "SHA256"], stage_rows) +
        key_values(execution, ["final_source_path", "final_source_sha256"]) + "</details>" +
        strategy_html + "<details><summary>Medaka model selection</summary>" + model_html + "</details></div></section>" +
        validation_section +
        "<section id='aware'><h2>AWARE decision explorer</h2><div class='panel'>"
        "<p class='note'>Target election, per-platform read actions and conservative confidence are separate outputs. "
        "A low confidence label does not override a platform's passthrough.</p>" +
        decision_notice + key_values(dec, keys) + "</div></section>"
        "<section id='removed'><h2>Removed-read origins</h2><div class='panel'>" + audit_html +
        f"<h3>{platform_label} filtering</h3>" + ont_html +
        "<p class='note'>Kraken classifications are taxonomic assignments, not confirmed biological origin. "
        "Paired reads are counted once per normalized pair identifier. No read sequences or IDs are exported.</p>"
        "</div></section>"
        "<section id='provenance'><h2>Evidence and exports</h2><div class='panel'>" +
        compute_html + "<h3>Discovered evidence</h3>" + table(["Evidence", "Source"], provenance_rows) +
        file_links +
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
    p.add_argument("--trajectory-dir", type=Path)
    p.add_argument("--validation-dir", type=Path)
    audit_group = p.add_mutually_exclusive_group()
    audit_group.add_argument("--audit-from", type=Path, help="Restore saved filtering audit records for a new report edition")
    audit_group.add_argument("--audit-reads", action="store_true",
                   help="Disk-backed comparison of merged/cleaned FASTQ IDs and saved Kraken assignments")
    args = p.parse_args()
    if not args.isolate or args.isolate in (".", "..") or "/" in args.isolate or "\\" in args.isolate:
        p.error("Invalid isolate identifier")
    iso = args.base_dir / args.isolate
    if not iso.is_dir():
        p.error("Isolate directory does not exist: " + str(iso))
    output = args.outdir or iso / "reports/hybas_research"
    result = generate(iso, output, args.audit_reads, args.trajectory_dir, args.validation_dir, args.audit_from)
    print("Report:", output / "hybas_report.html")
    print("Read audit:", result["status"])


if __name__ == "__main__":
    main()

"""Standalone HybAs report regression; no real Kraken/assembly required."""
import gzip
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/hybas_research_report.py"


def fq(path, names):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as f:
        for name in names:
            f.write(f"@{name}\nACGT\n+\nIIII\n")


def fixture(tmp_path, cleaned=("target",), decision="whitelist_filter"):
    iso = tmp_path / "isolate001"
    merged = iso / "work/inputs_merged"
    out = iso / "work/inputs_cleaned"
    for mate in ("R1", "R2"):
        fq(merged / f"illumina_{mate}.merged.fq.gz", ["target", "human", "other", "unknown"])
        fq(out / f"illumina_{mate}.clean.fq.gz", cleaned)
    kr = iso / "reports/kraken2"
    kr.mkdir(parents=True)
    (kr / "illumina.kraken2.txt").write_text(
        "C\ttarget\t1763\t4\tA\nC\thuman\t9606\t4\tA\n"
        "C\tother\t562\t4\tA\nU\tunknown\t0\t4\tA\n")
    (kr / "illumina.report.txt").write_text(
        "100.00\t4\t0\tR\t1\troot\n"
        " 75.00\t3\t0\tD\t2\t  Bacteria\n"
        " 25.00\t1\t1\tG\t1763\t    Mycobacterium\n"
        " 25.00\t1\t1\tS\t562\t    Escherichia coli\n"
        " 25.00\t1\t0\tD\t2759\t  Eukaryota\n"
        " 25.00\t1\t1\tS\t9606\t    Homo sapiens\n"
        " 25.00\t1\t1\tU\t0\t  unclassified\n")
    meta = iso / "metadata"
    meta.mkdir()
    (meta / "kraken_decision.tsv").write_text(
        "field\tvalue\nselected_target_taxid\t1763\n"
        "illumina_kraken_decision\t" + decision + "\n"
        "decision_confidence\tlow\nawareness_state\tfail_no_target\n")
    (meta / "run_summary.json").write_text(json.dumps({"run_mode": "illumina_only", "assembly_bp": 5000000}))
    (meta / "polish_execution.json").write_text(json.dumps({"executed_stages": ["polypolish"],
                                                              "final_source_stage": "polypolish"}))
    return iso


def run(tmp, *args):
    return subprocess.run([sys.executable, str(SCRIPT), "--base-dir", str(tmp),
                           "--isolate", "isolate001", *args],
                          capture_output=True, text=True, check=True)


def test_taxonomic_removed_read_reconciliation(tmp_path):
    iso = fixture(tmp_path)
    run(tmp_path, "--audit-reads")
    out = iso / "reports/hybas_research"
    result = json.loads((out / "aware_filter_audit.json").read_text())
    assert result["status"] == "completed"
    assert result["integrity_ok"]
    assert result["action_consistent"]
    assert result["totals"]["input_pairs"] == 4
    assert result["totals"]["removed_pairs"] == 3
    by_cat = {x["category"]: x for x in result["by_category"]}
    assert by_cat["Human-assigned"]["removed_pairs"] == 1
    assert by_cat["Other bacteria"]["removed_pairs"] == 1
    assert by_cat["Unclassified"]["removed_pairs"] == 1
    assert by_cat["Selected bacterial lineage"]["retained_pairs"] == 1
    page = (out / "hybas_report.html").read_text()
    assert "AWARE decision explorer" in page and "Human-assigned" in page
    assert "fail_no_target" in page and "polypolish" in page
    assert (out / "aware_removed_taxa.tsv").exists()


def test_passthrough_no_false_removals(tmp_path):
    fixture(tmp_path, cleaned=("target", "human", "other", "unknown"), decision="passthrough")
    run(tmp_path, "--audit-reads")
    result = json.loads((tmp_path / "isolate001/reports/hybas_research/aware_filter_audit.json").read_text())
    assert result["totals"]["removed_pairs"] == 0
    assert result["action_consistent"]


def test_partial_mate_retention_is_flagged(tmp_path):
    iso = fixture(tmp_path)
    fq(iso / "work/inputs_cleaned/illumina_R2.clean.fq.gz", [])
    run(tmp_path, "--audit-reads")
    result = json.loads((iso / "reports/hybas_research/aware_filter_audit.json").read_text())
    assert not result["integrity_ok"]
    assert result["totals"]["partial_pairs"] == 1


def test_report_without_read_audit(tmp_path):
    iso = fixture(tmp_path)
    run(tmp_path)
    out = iso / "reports/hybas_research"
    result = json.loads((out / "aware_filter_audit.json").read_text())
    assert result["status"] == "not_requested"
    assert (out / "hybas_report.html").is_file()


def test_ont_taxonomic_removed_reads(tmp_path):
    iso = fixture(tmp_path)
    before = iso / "work/inputs_merged/ont_merged.fastq.gz"
    after = iso / "work/inputs_cleaned/ont_clean.fastq.gz"
    fq(before, ["target", "human", "other"])
    fq(after, ["target"])
    kr = iso / "reports/kraken2"
    (kr / "nanopore.kraken2.txt").write_text(
        "C\ttarget\t1763\t4\tA\nC\thuman\t9606\t4\tA\nC\tother\t562\t4\tA\n")
    (kr / "nanopore.report.txt").write_text((kr / "illumina.report.txt").read_text())
    run(tmp_path, "--audit-reads")
    result = json.loads((iso / "reports/hybas_research/aware_ont_filter_audit.json").read_text())
    assert result["integrity_ok"]
    assert result["totals"]["removed_reads"] == 2
    categories = {row["category"]: row["removed_reads"] for row in result["by_category"]}
    assert categories["Human-assigned"] == 1
    assert categories["Other bacteria"] == 1
    assert (iso / "reports/hybas_research/aware_ont_removed_taxa.tsv").exists()


def test_historical_work_layout_and_stage_qc(tmp_path):
    iso = fixture(tmp_path)
    work = iso / "work"
    for mate in ("R1", "R2"):
        for suffix, old_dir in (("merged", "inputs_merged"), ("clean", "inputs_cleaned")):
            old = work / old_dir / f"illumina_{mate}.{suffix}.fq.gz"
            old.rename(work / old.name)
    (iso / "metadata/kraken_decision.tsv").unlink()
    quast = iso / "reports/quast"
    quast.mkdir()
    (quast / "report.txt").write_text(
        "Assembly\tassembly\n# contigs\t3\nTotal length\t4567890\nN50\t250000\nGC (%)\t62.5\n")
    stage = iso / "reports/FINAL_AUDIT_RESULTS/QC_ASSEMBLY/qc_steps/04_Polypolish"
    stage.mkdir(parents=True)
    (stage / "report.txt").write_text(
        "Assembly\tassembly\n# contigs\t2\nTotal length\t4567890\nN50\t450000\nGC (%)\t62.5\n")
    run(tmp_path, "--audit-reads")
    output = iso / "reports/hybas_research"
    result = json.loads((output / "aware_filter_audit.json").read_text())
    assert result["status"] == "completed"
    assert result["totals"]["removed_pairs"] == 3
    page = (output / "hybas_report.html").read_text()
    assert "4567890" in page and "450000" in page
    assert "Historical run: AWARE decision metadata was not found" in page
    assert "stage-graph" in page

###############################################################################
# Snakefile – HybAs 9.0.0-AWARE (native Snakemake reliability update)
#
# IEEE-grade design goals
# -----------------------
# 1. Batch mode across isolates via samples.tsv
# 2. Early config and input validation
# 3. Optional module toggles
# 4. Explicit per-isolate run mode declaration
# 5. Deterministic, stage-resolved polishing outputs
# 6. Structured provenance capture
# 7. Full environment manifest capture
# 8. Tool version and system provenance capture
# 9. Checksums for major inputs/outputs
# 10. Benchmarks separated from analytical outputs
# 11. Integrated summary TSV/JSON per isolate
# 12. Snakemake report-friendly outputs
# 13. Centralized databases under {base_dir}/db
# 14. Manual database update rules (never auto-download during normal runs)
# 15. Kraken decision layer with explicit off/manual/auto modes
# 16. Hardened Polypolish SAM generation
#
# Supported run modes
# -------------------
# - hybrid
# - ont_only
# - illumina_only
#
# Supported polish modes
# ----------------------
# - none
# - racon
# - medaka
# - tripolish
#
# "tripolish" corresponds to the named three-stage cascade:
#   Racon (configured rounds; manuscript default = 2)
#   -> Medaka
#   -> Polypolish
###############################################################################

import os
import csv
import json
import glob
import shlex
import hashlib
import platform
import re
import sys
from pathlib import Path

SOURCE_ROOT = os.path.abspath(workflow.basedir)
sys.path.insert(0, os.path.join(SOURCE_ROOT, "workflow", "lib"))
from workflow_support import (VERSION, normalize_taxid, manual_taxid, validate_isolate_names,
    pair_illumina_files, validate_input_files, polish_plan, busco_db_ready as valid_busco_db,
    file_identity, require_fasta, write_json, quast_metric as read_quast_metric)
from datetime import datetime
from computational import install_recorder, render_report
from settings import SettingsResolution, resource_value
from input_manifest import load_manifest
from library_inputs import execution_samples
from snakemake.io import report, directory, touch

def as_bool(x):
    if isinstance(x, bool):
        return x
    if isinstance(x, str):
        v = x.strip().lower()
        if v in {"true", "1", "yes", "y", "on"}:
            return True
        if v in {"false", "0", "no", "n", "off"}:
            return False
    raise ValueError(f"Invalid boolean {x!r}; use true or false")

configfile: "config.yaml"
report: "workflow/report/workflow.rst"

# Resolve before constructing rule parameters. Child jobs reuse the controller's
# frozen record; parsing/dry runs only preview and never write a settings file.
REQUESTED_CONFIG = json.dumps(config, sort_keys=True, default=str)
SETTINGS = SettingsResolution(config,
    cores=resource_value(workflow.global_resources, "_cores"),
    mem_mb=resource_value(workflow.global_resources, "mem_mb"))
config.clear()
config.update(SETTINGS.config)
workflow.global_resources["mem_mb"] = SETTINGS.memory_budget
workflow.global_resources["_cores"] = SETTINGS.cpu_budget
if workflow.is_main_process:
    print(SETTINGS.summary())

###############################################################################
# SECTION 1 — GLOBAL CONFIG
###############################################################################

BASE_DIR                 = os.path.abspath(os.path.expanduser(config.get("base_dir", ".")))
SAMPLE_SHEET             = config.get("sample_sheet", "")
INPUT_MANIFEST = load_manifest(os.path.expanduser(config["input_manifest"])) if config.get("input_manifest") else None
LIBRARY_SAMPLES = execution_samples(INPUT_MANIFEST) if INPUT_MANIFEST else {}
if INPUT_MANIFEST and INPUT_MANIFEST["source"]["sha256"] != SETTINGS.record["request_context"]["input_manifest_sha256"]:
    raise ValueError("Input manifest changed during configuration; retry with a stable manifest")
LIBRARY_ADAPTER = os.path.join(SOURCE_ROOT, "workflow", "lib", "library_inputs.py")
THREADS                  = int(config.get("threads", 8))
KEEP_PERCENT             = int(config.get("keep_percent", 95))
RACON_ROUNDS             = int(config.get("racon_rounds", 2))
DEFAULT_BUSCO_LINEAGE    = str(config.get("busco_lineage", "bacteria_odb10"))
DEFAULT_MEDAKA_MODEL     = str(config.get("medaka_model", ""))
DEFAULT_TARGET_TAXID     = manual_taxid(config.get("target_taxid"))

RUN_KRAKEN               = as_bool(config.get("run_kraken", True))
RUN_BUSCO                = as_bool(config.get("run_busco", True))
RUN_PROKKA               = as_bool(config.get("run_prokka", True))
RUN_MULTIQC              = as_bool(config.get("run_multiqc", True))
RUN_POLYPOLISH           = as_bool(config.get("run_polypolish", True))

SETUP_ONLY               = as_bool(config.get("setup_only", False))
ADMIN_MODE               = as_bool(config.get("admin_mode", False))
CI_MODE                  = as_bool(config.get("ci_mode", False))

ALLOW_ONT_ONLY           = as_bool(config.get("allow_ont_only", True))
ALLOW_ILLUMINA_ONLY      = as_bool(config.get("allow_illumina_only", True))
REQUIRE_TAXID_FOR_KRAKEN = as_bool(config.get("require_target_taxid_for_kraken", False))

POLISH_MODE              = str(config.get("polish_mode", "tripolish")).strip().lower()
WORKFLOW_VERSION         = VERSION
CONFIG_VERSION_LABEL     = str(config.get("workflow_version", VERSION))

DB_ROOT                  = os.path.join(BASE_DIR, config.get("db_root", "db"))
KRAKEN_DB_NAME           = str(config.get("kraken_db_name", "standard"))
KRAKEN_DB_DIR            = os.path.join(DB_ROOT, KRAKEN_DB_NAME)
BUSCO_DB_DIR             = os.path.join(DB_ROOT, "busco")
BUSCO_DOWNLOAD_LINEAGES  = list(config.get("busco_download_lineages", [DEFAULT_BUSCO_LINEAGE]))

# Kraken decision layer
KRAKEN_DECISION_MODE                 = str(config.get("kraken_decision_mode", "manual")).strip().lower()
KRAKEN_AUTO_PASSTHROUGH_MIN_PCT      = float(config.get("kraken_auto_passthrough_min_percent", 95.0))
KRAKEN_AUTO_USE_BOTH_PLATFORMS       = as_bool(config.get("kraken_auto_use_both_platforms", True))
KRAKEN_AUTO_PREFER_SPECIES_MIN_PCT   = float(config.get("kraken_auto_prefer_species_min_percent", 90.0))
KRAKEN_INCLUDE_CHILDREN              = as_bool(config.get("kraken_include_children", True))

# Awareness layer
KRAKEN_AWARE_DELTA_MIN_PCT           = float(config.get("kraken_aware_delta_min_percent", 20.0))
KRAKEN_AWARE_MODERATE_DELTA_MIN_PCT  = float(config.get("kraken_aware_moderate_delta_min_percent", 10.0))
KRAKEN_AWARE_HIGH_TARGET_MIN_PCT     = float(config.get("kraken_aware_high_target_min_percent", 95.0))
KRAKEN_AWARE_MODERATE_TARGET_MIN_PCT = float(config.get("kraken_aware_moderate_target_min_percent", 80.0))

VALID_POLISH_MODES = {"none", "racon", "medaka", "tripolish"}
VALID_KRAKEN_DECISION_MODES = {"off", "manual", "auto", "aware"}

# Effective Kraken profiling policy:
# - off    -> Kraken does not run
# - manual -> Kraken runs
# - auto   -> Kraken runs
# - aware  -> Kraken runs
KRAKEN_PROFILE_ENABLED = (KRAKEN_DECISION_MODE in {"manual", "auto", "aware"})


if POLISH_MODE not in VALID_POLISH_MODES:
    raise ValueError(f"Invalid polish_mode={POLISH_MODE!r}. Allowed: {sorted(VALID_POLISH_MODES)}")

if KRAKEN_DECISION_MODE not in VALID_KRAKEN_DECISION_MODES:
    raise ValueError(
        f"Invalid kraken_decision_mode={KRAKEN_DECISION_MODE!r}. "
        f"Allowed: {sorted(VALID_KRAKEN_DECISION_MODES)}"
    )

if THREADS < 1:
    raise ValueError("threads must be >= 1")
if not (1 <= KEEP_PERCENT <= 100):
    raise ValueError("keep_percent must be between 1 and 100")
if RACON_ROUNDS < 0:
    raise ValueError("racon_rounds must be >= 0")
if not (0.0 <= KRAKEN_AUTO_PASSTHROUGH_MIN_PCT <= 100.0):
    raise ValueError("kraken_auto_passthrough_min_percent must be between 0 and 100")

if not (0.0 <= KRAKEN_AUTO_PREFER_SPECIES_MIN_PCT <= 100.0):
    raise ValueError("kraken_auto_prefer_species_min_percent must be between 0 and 100")

if not (0.0 <= KRAKEN_AWARE_DELTA_MIN_PCT <= 100.0):
    raise ValueError("kraken_aware_delta_min_percent must be between 0 and 100")

if not (0.0 <= KRAKEN_AWARE_MODERATE_DELTA_MIN_PCT <= 100.0):
    raise ValueError("kraken_aware_moderate_delta_min_percent must be between 0 and 100")

if not (0.0 <= KRAKEN_AWARE_HIGH_TARGET_MIN_PCT <= 100.0):
    raise ValueError("kraken_aware_high_target_min_percent must be between 0 and 100")

if not (0.0 <= KRAKEN_AWARE_MODERATE_TARGET_MIN_PCT <= 100.0):
    raise ValueError("kraken_aware_moderate_target_min_percent must be between 0 and 100")

if KRAKEN_DECISION_MODE == "off" and RUN_KRAKEN:
    print("[WARN] kraken_decision_mode='off' overrides run_kraken=true; Kraken profiling will be skipped.")

if KRAKEN_DECISION_MODE in {"manual", "auto", "aware"} and not RUN_KRAKEN:
    raise ValueError(
        "kraken_decision_mode is 'manual' or 'auto' but run_kraken=false. "
        "Set run_kraken=true or switch kraken_decision_mode='off'."
    )
    
if not SETUP_ONLY and not ADMIN_MODE and not SAMPLE_SHEET and not LIBRARY_SAMPLES and "isolate" not in config:
    raise ValueError("Provide input_manifest, sample_sheet or legacy isolate in config.yaml")

###############################################################################
# SECTION 2 — DATABASE HELPERS
###############################################################################

def kraken_db_ready():
    expected = ["hash.k2d", "opts.k2d", "taxo.k2d"]
    return all(os.path.isfile(os.path.join(KRAKEN_DB_DIR, f)) and os.path.getsize(os.path.join(KRAKEN_DB_DIR, f)) > 0 for f in expected)

def busco_lineage_dir(lineage):
    return os.path.join(BUSCO_DB_DIR, lineage)

def busco_db_ready(lineage):
    return valid_busco_db(busco_lineage_dir(lineage))

###############################################################################
# SECTION 3 — SAMPLE SHEET LOADING
###############################################################################

def safe_float(x, default=0.0):
    try:
        return float(str(x).strip())
    except Exception:
        return default

def read_sample_sheet(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "isolate" not in reader.fieldnames:
            raise ValueError("samples.tsv must contain column 'isolate'")
        for row in reader:
            iso = str(row.get("isolate", "")).strip()
            if not iso:
                continue
            rows.append({
                "isolate": iso,
                "target_taxid": manual_taxid(row.get("target_taxid")),
                "busco_lineage": str(row.get("busco_lineage", "") or "").strip(),
                "medaka_model": str(row.get("medaka_model", "") or "").strip(),
            })
    return rows

if (SETUP_ONLY or ADMIN_MODE) and not LIBRARY_SAMPLES and (not SAMPLE_SHEET or not os.path.exists(SAMPLE_SHEET)):
    SAMPLE_ROWS = []
elif SAMPLE_SHEET:
    if not os.path.exists(SAMPLE_SHEET):
        raise ValueError(f"sample_sheet not found: {SAMPLE_SHEET}")
    SAMPLE_ROWS = read_sample_sheet(SAMPLE_SHEET)
elif LIBRARY_SAMPLES:
    SAMPLE_ROWS = [{"isolate": name, "target_taxid": "", "busco_lineage": "", "medaka_model": ""}
                   for name in LIBRARY_SAMPLES]
else:
    SAMPLE_ROWS = [{
        "isolate": str(config["isolate"]).strip(),
        "target_taxid": normalize_taxid(config.get("target_taxid", "")),
        "busco_lineage": str(config.get("busco_lineage", "") or "").strip(),
        "medaka_model": str(config.get("medaka_model", "") or "").strip(),
    }]

if not SAMPLE_ROWS and not (SETUP_ONLY or ADMIN_MODE):
    raise ValueError("No isolates found in sample sheet/config")

SAMPLES = [r["isolate"] for r in SAMPLE_ROWS]
validate_isolate_names(SAMPLES)
SAMPLE_INFO = {r["isolate"]: r for r in SAMPLE_ROWS}
if LIBRARY_SAMPLES and set(SAMPLES) != set(LIBRARY_SAMPLES):
    raise ValueError("sample_sheet and input_manifest must contain exactly the same sample IDs; no samples will be omitted")
LIBRARY_PATTERN = "(?:" + "|".join(re.escape(i) for i in LIBRARY_SAMPLES) + ")" if LIBRARY_SAMPLES else r"(?!)"
LEGACY_PATTERN = "(?:" + "|".join(re.escape(i) for i in SAMPLES if i not in LIBRARY_SAMPLES) + ")" if set(SAMPLES) - set(LIBRARY_SAMPLES) else r"(?!)"

###############################################################################
# SECTION 4 — PATH HELPERS
###############################################################################

def iso_dir(i): return os.path.join(BASE_DIR, i)

def raw_illumina_dir(i): return os.path.join(iso_dir(i), "raw", "illumina")
def raw_nanopore_dir(i): return os.path.join(iso_dir(i), "raw", "nanopore")

def work_dir(i): return os.path.join(iso_dir(i), "work")
def work_inputs_merged_dir(i): return os.path.join(work_dir(i), "inputs_merged")
def work_inputs_cleaned_dir(i): return os.path.join(work_dir(i), "inputs_cleaned")
def work_inputs_filtered_dir(i): return os.path.join(work_dir(i), "inputs_filtered")
def work_assembly_dir(i): return os.path.join(work_dir(i), "assembly")
def work_polishing_dir(i): return os.path.join(work_dir(i), "polishing")
def work_racon_dir(i): return os.path.join(work_polishing_dir(i), "racon")
def work_medaka_dir(i): return os.path.join(work_polishing_dir(i), "medaka")
def work_polypolish_dir(i): return os.path.join(work_polishing_dir(i), "polypolish")
def work_mapping_dir(i): return os.path.join(work_dir(i), "mapping")
def work_final_dir(i): return os.path.join(work_dir(i), "final")
def trimmed_dir(i): return os.path.join(iso_dir(i), "trimmed")

def reports_dir(i): return os.path.join(iso_dir(i), "reports")
def kraken_dir(i): return os.path.join(reports_dir(i), "kraken2")
def fastqc_pre_dir(i): return os.path.join(reports_dir(i), "fastqc_pre")
def fastp_dir(i): return os.path.join(reports_dir(i), "fastp")
def nanoplot_raw_dir(i): return os.path.join(reports_dir(i), "nanoplot_raw")
def nanoplot_filt_dir(i): return os.path.join(reports_dir(i), "nanoplot_filt")
def quast_dir(i): return os.path.join(reports_dir(i), "quast")
def busco_dir(i): return os.path.join(reports_dir(i), "busco")
def busco_run_dir(i): return os.path.join(busco_dir(i), f"busco_{i}")
def multiqc_dir(i): return os.path.join(reports_dir(i), "multiqc")
def seqkit_dir(i): return os.path.join(reports_dir(i), "seqkit")
def coverage_dir(i): return os.path.join(reports_dir(i), "coverage")

def annotation_dir(i): return os.path.join(iso_dir(i), "annotation")
def logs_dir(i): return os.path.join(iso_dir(i), "logs")
def metadata_dir(i): return os.path.join(iso_dir(i), "metadata")
def benchmarks_dir(i): return os.path.join(iso_dir(i), "benchmarks")

# merged / cleaned / filtered inputs
def merged_r1(i): return os.path.join(work_inputs_merged_dir(i), "illumina_R1.merged.fq.gz")
def merged_r2(i): return os.path.join(work_inputs_merged_dir(i), "illumina_R2.merged.fq.gz")
def merged_ont(i): return os.path.join(work_inputs_merged_dir(i), "ont_merged.fastq.gz")

def clean_r1(i): return os.path.join(work_inputs_cleaned_dir(i), "illumina_R1.clean.fq.gz")
def clean_r2(i): return os.path.join(work_inputs_cleaned_dir(i), "illumina_R2.clean.fq.gz")
def clean_ont_path(i): return os.path.join(work_inputs_cleaned_dir(i), "ont_clean.fastq.gz")

def trim_r1(i): return os.path.join(trimmed_dir(i), "R1.trimmed.fq.gz")
def trim_r2(i): return os.path.join(trimmed_dir(i), "R2.trimmed.fq.gz")
def filt_ont(i): return os.path.join(work_inputs_filtered_dir(i), "ont_filtered.fastq.gz")

# assembly / polishing / mapping / final
def unicycler_outdir(i): return os.path.join(work_assembly_dir(i), "unicycler")
def unicycler_asm(i): return os.path.join(unicycler_outdir(i), "assembly.fasta")
def stage_unicycler(i): return os.path.join(work_assembly_dir(i), "stage_unicycler.fasta")
def library_audit(i): return os.path.join(work_inputs_merged_dir(i), "library_identity")
def library_trimmed(i): return os.path.join(trimmed_dir(i), "libraries")
def spades_outdir(i): return os.path.join(work_assembly_dir(i), "spades_libraries")
def stage_spades(i): return os.path.join(work_assembly_dir(i), "stage_spades.fasta")
def base_stage(i): return "spades" if i in LIBRARY_SAMPLES else "unicycler"

def stage_racon0(i): return os.path.join(work_racon_dir(i), "stage_racon0.fasta")
def stage_racon(i, n): return os.path.join(work_racon_dir(i), f"stage_racon{n}.fasta")
def racon_paf(i, n): return os.path.join(work_racon_dir(i), f"ont_vs_racon{n}.paf")

def stage_medaka(i): return os.path.join(work_medaka_dir(i), "consensus.fasta")

def bwa_ref(i): return os.path.join(work_mapping_dir(i), "assembly.medaka.fa")
def bwa_fai(i): return os.path.join(work_mapping_dir(i), "assembly.medaka.fa.fai")
def bwa_idx(i, ext): return os.path.join(work_mapping_dir(i), f"assembly.medaka.fa.{ext}")
def coord_bam(i): return os.path.join(work_mapping_dir(i), "illumina.bam")
def coord_bai(i): return os.path.join(work_mapping_dir(i), "illumina.bam.bai")
def namesort_bam(i): return os.path.join(work_mapping_dir(i), "illumina.namesort.bam")

def stage_polypolish(i): return os.path.join(work_polypolish_dir(i), "stage_polypolish.fasta")
def final_asm(i): return os.path.join(work_final_dir(i), "assembly.final.fasta")

# reports
def quast_txt(i): return os.path.join(quast_dir(i), "report.txt")
def quast_tsv(i): return os.path.join(quast_dir(i), "report.tsv")
def busco_txt(i): return os.path.join(busco_run_dir(i), "short_summary.txt")
def busco_html(i): return os.path.join(busco_run_dir(i), "short_summary.html")
def multiqc_html(i): return os.path.join(multiqc_dir(i), "multiqc_report.html")

# annotation
def prokka_gff(i): return os.path.join(annotation_dir(i), f"{i}.gff")
def prokka_faa(i): return os.path.join(annotation_dir(i), f"{i}.faa")

# metadata
def run_mode_file(i): return os.path.join(metadata_dir(i), "run_mode.txt")
def run_summary_tsv(i): return os.path.join(metadata_dir(i), "run_summary.tsv")
def run_summary_json(i): return os.path.join(metadata_dir(i), "run_summary.json")
def qc_assessment_tsv(i): return os.path.join(metadata_dir(i), "qc_assessment.tsv")
def multiqc_custom_tsv(i): return os.path.join(metadata_dir(i), "multiqc_custom_content.tsv")
def system_info_tsv(i): return os.path.join(metadata_dir(i), "system_info.tsv")
def workflow_provenance_tsv(i): return os.path.join(metadata_dir(i), "workflow_provenance.tsv")
def tool_versions_tsv(i): return os.path.join(metadata_dir(i), "tool_versions.tsv")
def checksums_tsv(i): return os.path.join(metadata_dir(i), "checksums.tsv")
def skip_reasons_tsv(i): return os.path.join(metadata_dir(i), "skip_reasons.tsv")
def kraken_decision_tsv(i): return os.path.join(metadata_dir(i), "kraken_decision.tsv")


# seqkit / coverage
def seqkit_illumina(i): return os.path.join(seqkit_dir(i), "illumina_pre.tsv")
def seqkit_ont(i): return os.path.join(seqkit_dir(i), "ont.tsv")
def seqkit_assembly(i): return os.path.join(seqkit_dir(i), "assembly.tsv")
def coverage_note(i): return os.path.join(coverage_dir(i), "ONT_coverage.txt")
def coverage_checked(i): return os.path.join(coverage_dir(i), "ONT_coverage.checked")

# benchmarks
def bench(i, name): return os.path.join(benchmarks_dir(i), f"{name}.tsv")

# Kraken aware mode assembly type selector
def combine_platform_values(val_illumina, val_ont, use_both=True):
    vals = []
    try:
        if val_illumina is not None:
            vals.append(("illumina", float(val_illumina)))
    except Exception:
        pass
    try:
        if val_ont is not None:
            vals.append(("ont", float(val_ont)))
    except Exception:
        pass

    # keep only meaningful non-negative values from present platforms
    vals = [(name, v) for name, v in vals if v >= 0]

    if not vals:
        return 0.0, "no_platform_values"

    if len(vals) == 1:
        return vals[0][1], f"single_platform_{vals[0][0]}"

    if use_both:
        return min(v for _, v in vals), "both_platforms_min"
    return max(v for _, v in vals), "both_platforms_max"

def choose_platform_cleaning_decision(mode, selected_taxid, platform_pct, platform_name):
    selected_taxid = normalize_taxid(selected_taxid)

    if mode == "off":
        return ("passthrough", f"{platform_name}_mode_off")

    if mode == "manual":
        if not selected_taxid:
            return ("passthrough", f"{platform_name}_manual_mode_no_target_taxid")
        return ("whitelist_filter", f"{platform_name}_manual_mode_target_taxid_supplied")

    # auto / aware: decide per platform
    if not selected_taxid:
        return ("passthrough", f"{platform_name}_{mode}_mode_no_selected_target_taxid")

    try:
        pct = float(platform_pct)
    except Exception:
        pct = 0.0

    if pct >= KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:
        return (
            "passthrough",
            f"{platform_name}_{mode}_mode_selected_target_pct={pct:.2f}; "
            f"threshold={KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:.2f}"
        )

    return (
        "whitelist_filter",
        f"{platform_name}_{mode}_mode_selected_target_pct={pct:.2f}; "
        f"threshold={KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:.2f}"
    )


###############################################################################
# PROVENANCE HELPERS
###############################################################################

PROJECT_CONDARC = os.path.abspath(os.environ.get("CONDARC", os.path.join(SOURCE_ROOT, ".condarc")))

def provenance_dir():
    return os.path.join(BASE_DIR, "provenance")

def provenance_conda_dir():
    return os.path.join(provenance_dir(), "conda")

def provenance_db_dir():
    return os.path.join(provenance_dir(), "databases")

def provenance_runtime_dir():
    return os.path.join(provenance_dir(), "runtime")

def provenance_condarc_file():
    return os.path.join(provenance_conda_dir(), "condarc.used.yml")

def provenance_system_info():
    return os.path.join(provenance_runtime_dir(), "system_info.tsv")

def provenance_workflow_info():
    return os.path.join(provenance_runtime_dir(), "workflow_provenance.tsv")

def provenance_db_hashes():
    return os.path.join(provenance_db_dir(), "database_hashes.tsv")

def provenance_input_hashes():
    return os.path.join(provenance_runtime_dir(), "global_input_hashes.tsv")

def env_conda_list_global(envname):
    return os.path.join(provenance_conda_dir(), f"{envname}.conda_list.txt")

def env_export_global(envname):
    return os.path.join(provenance_conda_dir(), f"{envname}.env_export.yml")

def env_explicit_global(envname):
    return os.path.join(provenance_conda_dir(), f"{envname}.explicit.txt")

def env_pip_freeze_global(envname):
    return os.path.join(provenance_conda_dir(), f"{envname}.pip_freeze.txt")

def env_python_version_global(envname):
    return os.path.join(provenance_conda_dir(), f"{envname}.python_version.txt")

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


    
###############################################################################
# METADATA / PARSING HELPERS
###############################################################################

def safe_read_text(path):
    if not os.path.exists(path):
        return ""
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        return handle.read()

def read_tsv_rows(path):
    txt = safe_read_text(path).strip()
    if not txt:
        return []
    lines = txt.splitlines()
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        vals = line.split("\t")
        if len(vals) < len(header):
            vals += [""] * (len(header) - len(vals))
        rows.append(dict(zip(header, vals[:len(header)])))
    return rows

def read_key_value_tsv(path, key_col=0, val_col=1, has_header=True):
    txt = safe_read_text(path).strip()
    if not txt:
        return {}
    lines = txt.splitlines()
    if has_header and len(lines) > 1:
        lines = lines[1:]
    out = {}
    for line in lines:
        parts = line.split("\t")
        if len(parts) > max(key_col, val_col):
            out[parts[key_col].strip()] = parts[val_col].strip()
    return out

def parse_seqkit_rows(path):
    txt = safe_read_text(path).strip()
    if not txt:
        return []
    lines = txt.splitlines()
    if len(lines) < 2:
        return []
    header = lines[0].split()
    rows = []
    for line in lines[1:]:
        vals = line.split()
        if len(vals) < len(header):
            vals += [""] * (len(header) - len(vals))
        rows.append(dict(zip(header, vals[:len(header)])))
    return rows

def seqkit_sum_len(row):
    for key in ("sum_len", "sum_len(bp)", "sum_len_bp"):
        if key in row and row[key]:
            return str(row[key]).replace(",", "")
    vals = list(row.values())
    if len(vals) >= 5:
        return str(vals[4]).replace(",", "")
    return "NA"

def parse_seqkit_total_bp_first_row(path):
    rows = parse_seqkit_rows(path)
    if not rows:
        return "NA"
    return seqkit_sum_len(rows[0])

def parse_seqkit_total_bp_all_rows(path):
    rows = parse_seqkit_rows(path)
    total = 0
    seen = False
    for row in rows:
        v = seqkit_sum_len(row)
        if v != "NA":
            try:
                total += int(v)
                seen = True
            except Exception:
                pass
    return str(total) if seen else "NA"

def parse_ont_raw_filtered_bp(path):
    rows = parse_seqkit_rows(path)
    raw_bp = "NA"
    filt_bp = "NA"
    if len(rows) >= 1:
        raw_bp = seqkit_sum_len(rows[0])
    if len(rows) >= 2:
        filt_bp = seqkit_sum_len(rows[1])
    return raw_bp, filt_bp

def fasta_contig_count(path):
    if not os.path.exists(path):
        return "NA"
    n = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                n += 1
    return str(n)

def parse_quast_report_tsv(path):
    """
    QUAST report.tsv is metric-transposed:
    Assembly <tab> assembly.final
    # contigs <tab> 1
    Total length <tab> 4719737
    GC (%) <tab> 37.19
    N50 <tab> 4719737
    """
    rows = read_tsv_rows(path)
    data = {}
    for row in rows:
        vals = list(row.values())
        if len(vals) >= 2:
            metric = vals[0].strip()
            value = vals[1].strip()
            if metric:
                data[metric] = value
    return data

def parse_busco_short_summary(path):
    txt = safe_read_text(path)
    out = {
        "summary_line": "NA",
        "complete_percent": "NA",
        "single_copy_percent": "NA",
        "duplicated_percent": "NA",
        "fragmented_percent": "NA",
        "missing_percent": "NA",
        "n_buscos": "NA",
        "complete_count": "NA",
        "single_copy_count": "NA",
        "duplicated_count": "NA",
        "fragmented_count": "NA",
        "missing_count": "NA",
    }

    import re

    for line in txt.splitlines():
        s = line.strip()

        if s.startswith("C:") or s.startswith("Complete BUSCOs"):
            out["summary_line"] = s

        m = re.search(
            r"C:(?P<C>[0-9.]+)%\[S:(?P<S>[0-9.]+)%,D:(?P<D>[0-9.]+)%\],F:(?P<F>[0-9.]+)%,M:(?P<M>[0-9.]+)%,n:(?P<n>[0-9]+)",
            s
        )
        if m:
            out["complete_percent"] = m.group("C")
            out["single_copy_percent"] = m.group("S")
            out["duplicated_percent"] = m.group("D")
            out["fragmented_percent"] = m.group("F")
            out["missing_percent"] = m.group("M")
            out["n_buscos"] = m.group("n")

        m2 = re.search(
            r"(?P<Cc>[0-9]+)\s+Complete BUSCOs.*?(?P<Sc>[0-9]+)\s+Complete and single-copy BUSCOs.*?(?P<Dc>[0-9]+)\s+Complete and duplicated BUSCOs.*?(?P<Fc>[0-9]+)\s+Fragmented BUSCOs.*?(?P<Mc>[0-9]+)\s+Missing BUSCOs.*?(?P<nc>[0-9]+)\s+Total BUSCO groups searched",
            txt.replace("\n", " ")
        )
        if m2:
            out["complete_count"] = m2.group("Cc")
            out["single_copy_count"] = m2.group("Sc")
            out["duplicated_count"] = m2.group("Dc")
            out["fragmented_count"] = m2.group("Fc")
            out["missing_count"] = m2.group("Mc")
            out["n_buscos"] = m2.group("nc")

    return out

def parse_cov_tsv(path):
    return read_key_value_tsv(path, key_col=0, val_col=1, has_header=True)

def safe_fraction(num, den):
    try:
        num = float(num)
        den = float(den)
        if den > 0:
            return f"{num/den:.4f}"
    except Exception:
        pass
    return "NA"

def parse_conda_list_version(conda_list_path, package_names):
    if isinstance(package_names, str):
        package_names = [package_names]

    package_names = {str(x).strip().lower() for x in package_names if str(x).strip()}
    if not os.path.exists(conda_list_path):
        return "NA"

    with open(conda_list_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            pkg = parts[0].strip().lower()
            ver = parts[1].strip()

            if pkg in package_names:
                return ver

    return "NA"

def parse_kraken_report_rows(report_path):
    rows = []
    if not os.path.exists(report_path):
        return rows

    with open(report_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rows.append({
                "pct": safe_float(parts[0], 0.0),
                "clade_reads": parts[1].strip(),
                "direct_reads": parts[2].strip(),
                "rank": parts[3].strip(),
                "taxid": normalize_taxid(parts[4]),
                "name": parts[5].strip(),
            })
    return rows

def parse_kraken_report_target_fraction(report_path, target_taxid):
    target_taxid = normalize_taxid(target_taxid)
    if not target_taxid:
        return 0.0

    for row in parse_kraken_report_rows(report_path):
        if row["taxid"] == target_taxid:
            return row["pct"]
    return 0.0

def is_meaningful_kraken_candidate(row):
    rank = row["rank"]
    name = row["name"].strip().lower()

    if not row["taxid"] or row["taxid"] == "0":
        return False
    if name.startswith("unclassified"):
        return False
    if rank in {"U", "R", "R1", "R2", "K", "P", "C", "O", "F"}:
        return False

    # allow only genus/species level primary nodes for auto target discovery
    return rank in {"G", "S"}

def infer_target_from_kraken_report(report_path):
    rows = [r for r in parse_kraken_report_rows(report_path) if is_meaningful_kraken_candidate(r)]

    if not rows:
        return {
            "taxid": "",
            "name": "NA",
            "rank": "NA",
            "pct": 0.0,
            "source_reason": "no_meaningful_candidate"
        }

    best_species = None
    best_genus = None

    for row in rows:
        if row["rank"] == "S":
            if best_species is None or row["pct"] > best_species["pct"]:
                best_species = row
        elif row["rank"] == "G":
            if best_genus is None or row["pct"] > best_genus["pct"]:
                best_genus = row

    if best_species is not None and best_species["pct"] >= KRAKEN_AUTO_PREFER_SPECIES_MIN_PCT:
        chosen = best_species
        why = f"species_preferred_at_{best_species['pct']:.2f}_percent"
    elif best_genus is not None:
        chosen = best_genus
        why = f"genus_preferred_at_{best_genus['pct']:.2f}_percent"
    elif best_species is not None:
        chosen = best_species
        why = f"species_fallback_at_{best_species['pct']:.2f}_percent"
    else:
        return {
            "taxid": "",
            "name": "NA",
            "rank": "NA",
            "pct": 0.0,
            "source_reason": "no_candidate_after_selection"
        }

    return {
        "taxid": chosen["taxid"],
        "name": chosen["name"],
        "rank": chosen["rank"],
        "pct": chosen["pct"],
        "source_reason": why
    }

def choose_selected_kraken_target(mode, manual_target_taxid, illumina_report, ont_report):
    manual_target_taxid = normalize_taxid(manual_target_taxid)

    if mode == "off":
        return {
            "selected_taxid": "",
            "selected_name": "NA",
            "selected_rank": "NA",
            "selected_source": "off_mode",
            "illumina_selected_pct": 0.0,
            "ont_selected_pct": 0.0,
            "selected_pct_for_decision": 0.0,
            "inferred_taxid": "",
            "inferred_name": "NA",
            "inferred_rank": "NA",
            "inferred_pct": 0.0,
            "selection_reason": "kraken_mode_off"
        }

    if mode == "manual":
        if not manual_target_taxid:
            return {
                "selected_taxid": "",
                "selected_name": "NA",
                "selected_rank": "NA",
                "selected_source": "manual_empty",
                "illumina_selected_pct": 0.0,
                "ont_selected_pct": 0.0,
                "selected_pct_for_decision": 0.0,
                "inferred_taxid": "",
                "inferred_name": "NA",
                "inferred_rank": "NA",
                "inferred_pct": 0.0,
                "selection_reason": "manual_mode_no_target_taxid"
            }

        illumina_pct = parse_kraken_report_target_fraction(illumina_report, manual_target_taxid)
        ont_pct = parse_kraken_report_target_fraction(ont_report, manual_target_taxid)

        metric, metric_mode = combine_platform_values(
            illumina_pct if os.path.isfile(illumina_report) and os.path.getsize(illumina_report) > 0 else None,
            ont_pct if os.path.isfile(ont_report) and os.path.getsize(ont_report) > 0 else None,
            use_both=KRAKEN_AUTO_USE_BOTH_PLATFORMS,
        )

        return {
            "selected_taxid": manual_target_taxid,
            "selected_name": "MANUAL_TARGET",
            "selected_rank": "MANUAL",
            "selected_source": "manual_target_taxid",
            "illumina_selected_pct": illumina_pct,
            "ont_selected_pct": ont_pct,
            "selected_pct_for_decision": metric,
            "inferred_taxid": "",
            "inferred_name": "NA",
            "inferred_rank": "NA",
            "inferred_pct": 0.0,
            "selection_reason": "manual_mode_target_taxid_supplied"
        }

    # auto and aware share the same target-selection logic:
    # 1) manual override if TSV taxid exists
    # 2) otherwise infer from Kraken
    if mode in {"auto", "aware"}:
        if manual_target_taxid:
            illumina_pct = parse_kraken_report_target_fraction(illumina_report, manual_target_taxid)
            ont_pct = parse_kraken_report_target_fraction(ont_report, manual_target_taxid)
            metric, metric_mode = combine_platform_values(
                illumina_pct if os.path.isfile(illumina_report) and os.path.getsize(illumina_report) > 0 else None,
                ont_pct if os.path.isfile(ont_report) and os.path.getsize(ont_report) > 0 else None,
                use_both=KRAKEN_AUTO_USE_BOTH_PLATFORMS,
            )

            return {
                "selected_taxid": manual_target_taxid,
                "selected_name": "MANUAL_TARGET_IN_AUTO",
                "selected_rank": "MANUAL",
                "selected_source": "manual_override_in_auto" if mode == "auto" else "manual_override_in_aware",
                "illumina_selected_pct": illumina_pct,
                "ont_selected_pct": ont_pct,
                "selected_pct_for_decision": metric,
                "inferred_taxid": "",
                "inferred_name": "NA",
                "inferred_rank": "NA",
                "inferred_pct": 0.0,
                "selection_reason": "auto_mode_manual_override_used" if mode == "auto" else "aware_mode_manual_override_used"
            }

        illumina_inferred = infer_target_from_kraken_report(illumina_report)
        ont_inferred = infer_target_from_kraken_report(ont_report)

        chosen = illumina_inferred if illumina_inferred["pct"] >= ont_inferred["pct"] else ont_inferred
        selected_taxid = normalize_taxid(chosen["taxid"])

        if not selected_taxid:
            return {
                "selected_taxid": "",
                "selected_name": "NA",
                "selected_rank": "NA",
                "selected_source": "auto_inference_failed" if mode == "auto" else "aware_inference_failed",
                "illumina_selected_pct": 0.0,
                "ont_selected_pct": 0.0,
                "selected_pct_for_decision": 0.0,
                "inferred_taxid": "",
                "inferred_name": "NA",
                "inferred_rank": "NA",
                "inferred_pct": 0.0,
                "selection_reason": "auto_mode_no_inferred_target" if mode == "auto" else "aware_mode_no_inferred_target"
            }

        illumina_pct = parse_kraken_report_target_fraction(illumina_report, selected_taxid)
        ont_pct = parse_kraken_report_target_fraction(ont_report, selected_taxid)
        metric = min(illumina_pct, ont_pct) if KRAKEN_AUTO_USE_BOTH_PLATFORMS else max(illumina_pct, ont_pct)

        return {
            "selected_taxid": selected_taxid,
            "selected_name": chosen["name"],
            "selected_rank": chosen["rank"],
            "selected_source": "inferred_from_kraken",
            "illumina_selected_pct": illumina_pct,
            "ont_selected_pct": ont_pct,
            "selected_pct_for_decision": metric,
            "inferred_taxid": chosen["taxid"],
            "inferred_name": chosen["name"],
            "inferred_rank": chosen["rank"],
            "inferred_pct": chosen["pct"],
            "selection_reason": chosen["source_reason"]
        }

    return {
        "selected_taxid": "",
        "selected_name": "NA",
        "selected_rank": "NA",
        "selected_source": "unhandled_mode",
        "illumina_selected_pct": 0.0,
        "ont_selected_pct": 0.0,
        "selected_pct_for_decision": 0.0,
        "inferred_taxid": "",
        "inferred_name": "NA",
        "inferred_rank": "NA",
        "inferred_pct": 0.0,
        "selection_reason": f"unhandled_mode_{mode}"
    }

def choose_kraken_decision(mode, selected_taxid, selected_pct_for_decision):
    selected_taxid = normalize_taxid(selected_taxid)

    if mode == "off":
        return ("passthrough", "kraken_mode_off")

    if mode == "manual":
        if not selected_taxid:
            return ("passthrough", "manual_mode_no_target_taxid")
        return ("whitelist_filter", "manual_mode_target_taxid_supplied")

    # auto / aware base decision
    if not selected_taxid:
        return ("passthrough", f"{mode}_mode_no_selected_target_taxid")

    if selected_pct_for_decision >= KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:
        return (
            "passthrough",
            f"{mode}_mode_selected_target_pct={selected_pct_for_decision:.2f}; "
            f"threshold={KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:.2f}"
        )

    return (
        "whitelist_filter",
        f"{mode}_mode_selected_target_pct={selected_pct_for_decision:.2f}; "
        f"threshold={KRAKEN_AUTO_PASSTHROUGH_MIN_PCT:.2f}"
    )

def top_two_meaningful_candidates(report_path):
    rows = [r for r in parse_kraken_report_rows(report_path) if is_meaningful_kraken_candidate(r)]
    rows = sorted(rows, key=lambda x: x["pct"], reverse=True)
    top1 = rows[0] if len(rows) >= 1 else None
    top2 = rows[1] if len(rows) >= 2 else None
    return top1, top2

def compute_platform_agreement(illumina_report, ont_report):
    i = infer_target_from_kraken_report(illumina_report)
    o = infer_target_from_kraken_report(ont_report)

    if not i["taxid"] or not o["taxid"]:
        return 0.0, "no_inference"

    if i["taxid"] == o["taxid"]:
        return 1.0, "exact_taxid_match"

    if i["rank"] == "S" and o["rank"] == "S":
        i_name = i["name"].split()[0] if i["name"] else ""
        o_name = o["name"].split()[0] if o["name"] else ""
        if i_name and i_name == o_name:
            return 0.5, "same_genus_different_species"

    return 0.0, "discordant_targets"

def compute_taxonomic_confidence(illumina_report, ont_report, selected_taxid):
    selected_taxid = normalize_taxid(selected_taxid)
    if not selected_taxid:
        return {
            "runner_up_taxid": "",
            "runner_up_name": "NA",
            "runner_up_pct": 0.0,
            "delta_pct": 0.0,
            "platform_agreement_score": 0.0,
            "platform_agreement_reason": "no_selected_target",
            "decision_confidence": "low",
            "awareness_state": "fail_no_target"
        }

    top1_i, top2_i = top_two_meaningful_candidates(illumina_report)
    top1_o, top2_o = top_two_meaningful_candidates(ont_report)

    selected_pct_i = parse_kraken_report_target_fraction(illumina_report, selected_taxid)
    selected_pct_o = parse_kraken_report_target_fraction(ont_report, selected_taxid)
    selected_pct = min(selected_pct_i, selected_pct_o) if KRAKEN_AUTO_USE_BOTH_PLATFORMS else max(selected_pct_i, selected_pct_o)

    runner_candidates = []
    for x in [top1_i, top2_i, top1_o, top2_o]:
        if x and normalize_taxid(x["taxid"]) != selected_taxid:
            runner_candidates.append(x)

    runner_candidates = sorted(runner_candidates, key=lambda x: x["pct"], reverse=True)
    runner = runner_candidates[0] if runner_candidates else None
    runner_pct = runner["pct"] if runner else 0.0
    delta_pct = max(selected_pct - runner_pct, 0.0)

    agreement_score, agreement_reason = compute_platform_agreement(illumina_report, ont_report)

    if selected_pct >= KRAKEN_AWARE_HIGH_TARGET_MIN_PCT and delta_pct >= KRAKEN_AWARE_DELTA_MIN_PCT and agreement_score >= 0.5:
        confidence = "high"
        state = "clean_stable"
    elif selected_pct >= KRAKEN_AWARE_MODERATE_TARGET_MIN_PCT and delta_pct >= KRAKEN_AWARE_MODERATE_DELTA_MIN_PCT:
        confidence = "moderate"
        state = "clean_filterable"
    elif selected_pct > 0:
        confidence = "low"
        state = "taxonomically_ambiguous"
    else:
        confidence = "low"
        state = "fail_no_target"

    return {
        "runner_up_taxid": normalize_taxid(runner["taxid"]) if runner else "",
        "runner_up_name": runner["name"] if runner else "NA",
        "runner_up_pct": runner_pct,
        "delta_pct": delta_pct,
        "platform_agreement_score": agreement_score,
        "platform_agreement_reason": agreement_reason,
        "decision_confidence": confidence,
        "awareness_state": state
    }

def normalize_tool_version_string(s):
    s = str(s).strip()
    if not s:
        return "NA"
    first = s.splitlines()[0].strip()
    if not first:
        return "NA"
    bad_prefixes = (
        "/usr/bin/bash:",
        "WARNING:",
        "Copyright ",
    )
    for bp in bad_prefixes:
        if first.startswith(bp):
            return "NA"
    return first


###############################################################################
# SECTION 9 — FINAL TARGETS
###############################################################################


def planned_stages(i):
    return polish_plan(RUN_MODE.get(i, "setup_only"), POLISH_MODE, RACON_ROUNDS, RUN_POLYPOLISH)

def needs_polypolish(i):
    return "polypolish" in planned_stages(i)

def stage_path(i, stage):
    if stage == "unicycler": return stage_unicycler(i)
    if stage == "spades": return stage_spades(i)
    if stage == "medaka": return stage_medaka(i)
    if stage == "polypolish": return stage_polypolish(i)
    if re.fullmatch(r"racon[0-9]+", stage): return stage_racon(i, int(stage[5:]))
    raise ValueError(f"Unknown assembly stage: {stage}")

def polishing_source(i):
    stages = [x for x in planned_stages(i) if x != "polypolish"]
    return stage_path(i, stages[-1] if stages else base_stage(i))

def polish_status(i): return os.path.join(metadata_dir(i), "polypolish_status.json")
def polish_execution(i): return os.path.join(metadata_dir(i), "polish_execution.json")
def polishing_sam(i, mate): return os.path.join(work_mapping_dir(i), f"polypolish_R{mate}.sam")
def alignment_count(i): return os.path.join(metadata_dir(i), "illumina_alignment_count.txt")

def enabled_envs():
    names = ["assembly", "qc"]
    if LIBRARY_SAMPLES: names.append("library_assembly")
    if (SETUP_ONLY and POLISH_MODE == "tripolish" and RUN_POLYPOLISH) or any(needs_polypolish(i) for i in SAMPLES): names.append("polishing")
    if (SETUP_ONLY and POLISH_MODE in {"medaka", "tripolish"}) or any("medaka" in planned_stages(i) for i in SAMPLES): names.append("medaka")
    if RUN_PROKKA: names.append("annotation")
    return names

def environment_outputs():
    return [fn(name) for name in enabled_envs() for fn in
            (env_conda_list_global, env_export_global, env_explicit_global,
             env_pip_freeze_global, env_python_version_global)]

def active_env_list(name):
    return env_conda_list_global(name) if name in enabled_envs() else []

def rule_memory(name, fallback):
    value = int(config.get(name + "_mem_mb", fallback))
    if value < 1: raise ValueError(f"{name}_mem_mb must be positive")
    return value

ASSEMBLY_MEM_MB = rule_memory("assembly", 8192)
SPADES_MEMORY_GB = int(config.get("spades_memory_gb", max(1, ASSEMBLY_MEM_MB // 1024)))
if SPADES_MEMORY_GB < 1 or SPADES_MEMORY_GB * 1024 > ASSEMBLY_MEM_MB:
    raise ValueError("spades_memory_gb must be positive and fit within assembly_mem_mb")
MIN_FASTA_LENGTH = int(config.get("min_fasta_length", 1000))
if MIN_FASTA_LENGTH < 1: raise ValueError("min_fasta_length must be positive")
KRAKEN_MEM_MB = rule_memory("kraken", max(1024, int(sum(
    os.path.getsize(os.path.join(KRAKEN_DB_DIR, f)) for f in ("hash.k2d", "opts.k2d", "taxo.k2d")
    if os.path.isfile(os.path.join(KRAKEN_DB_DIR, f))) / 1048576 * 1.2) + 1024))
POLISH_MEM_MB = rule_memory("polishing", 4096)
DB_CORE = [os.path.join(KRAKEN_DB_DIR, f) for f in ("hash.k2d", "opts.k2d", "taxo.k2d")]
CONFIG_FILES = [os.path.abspath(x) for x in list(workflow.configfiles) if os.path.isfile(x)]
SAMPLE_FILES = ([SAMPLE_SHEET] if SAMPLE_SHEET and os.path.isfile(SAMPLE_SHEET) else []) + ([INPUT_MANIFEST["source"]["manifest_path"]] if INPUT_MANIFEST else [])
SUPPORT_FILE = os.path.join(SOURCE_ROOT, "workflow", "lib", "workflow_support.py")
SOURCE_FILES = [os.path.join(SOURCE_ROOT, "Hybas.smk"), SUPPORT_FILE,
                os.path.join(SOURCE_ROOT, "workflow", "lib", "settings.py"), LIBRARY_ADAPTER,
                os.path.join(SOURCE_ROOT, "workflow", "lib", "input_manifest.py")]
EFFECTIVE_CONFIG = json.dumps(config, sort_keys=True, default=str)
AWARE_SETTINGS = {k: v for k, v in globals().copy().items() if k.startswith("KRAKEN_") and isinstance(v, (str, bool, int, float))}
# Explicitly fingerprint only controller source so report-only edits do not rerun classification.
_controller_source = Path(SOURCE_FILES[0]).read_text()
_controller_functions = ["safe_float", "combine_platform_values", "choose_platform_cleaning_decision"]
_controller_blocks = [re.search(r"(?ms)^def " + name + r"\(.*?(?=^def |^rule |\Z)", _controller_source).group(0)
                      for name in _controller_functions]
_controller_blocks.append(_controller_source[
    _controller_source.index("def parse_kraken_report_rows("):
    _controller_source.index("def normalize_tool_version_string(")])
_controller_blocks.append(Path(SUPPORT_FILE).read_text().split("def manual_taxid", 1)[0])
AWARE_CODE_HASH = hashlib.sha256("\n".join(_controller_blocks).encode()).hexdigest()

rule all:
    input:
        lambda wc: environment_outputs(),
        provenance_condarc_file(),
        os.path.join(provenance_runtime_dir(), "computational_report.html"),
        provenance_system_info(),
        provenance_workflow_info(),
        provenance_db_hashes(),
        provenance_input_hashes(),


        expand(final_asm("{isolate}"), isolate=SAMPLES),
        expand(polish_execution("{isolate}"), isolate=SAMPLES),
        expand(run_mode_file("{isolate}"), isolate=SAMPLES),
        expand(run_summary_tsv("{isolate}"), isolate=SAMPLES),
        expand(run_summary_json("{isolate}"), isolate=SAMPLES),
        expand(system_info_tsv("{isolate}"), isolate=SAMPLES),
        expand(workflow_provenance_tsv("{isolate}"), isolate=SAMPLES),
        expand(tool_versions_tsv("{isolate}"), isolate=SAMPLES),
        expand(checksums_tsv("{isolate}"), isolate=SAMPLES),
        expand(skip_reasons_tsv("{isolate}"), isolate=SAMPLES),
        expand(kraken_decision_tsv("{isolate}"), isolate=SAMPLES),
        expand(os.path.join(metadata_dir("{isolate}"), "qc_assessment.tsv"), isolate=SAMPLES),

        expand(fastqc_pre_dir("{isolate}"), isolate=SAMPLES),
        expand(nanoplot_raw_dir("{isolate}"), isolate=SAMPLES),
        expand(nanoplot_filt_dir("{isolate}"), isolate=SAMPLES),

        expand(os.path.join(kraken_dir("{isolate}"), "illumina.kraken2.txt"), isolate=SAMPLES),
        expand(os.path.join(kraken_dir("{isolate}"), "illumina.report.txt"), isolate=SAMPLES),
        expand(os.path.join(kraken_dir("{isolate}"), "nanopore.kraken2.txt"), isolate=SAMPLES),
        expand(os.path.join(kraken_dir("{isolate}"), "nanopore.report.txt"), isolate=SAMPLES),

        expand(os.path.join(fastp_dir("{isolate}"), "fastp.html"), isolate=SAMPLES),
        expand(os.path.join(fastp_dir("{isolate}"), "fastp.json"), isolate=SAMPLES),

        expand(seqkit_illumina("{isolate}"), isolate=SAMPLES),
        expand(seqkit_ont("{isolate}"), isolate=SAMPLES),
        expand(seqkit_assembly("{isolate}"), isolate=SAMPLES),
        expand(coverage_note("{isolate}"), isolate=SAMPLES),
        expand(coverage_checked("{isolate}"), isolate=SAMPLES),

        expand(quast_txt("{isolate}"), isolate=SAMPLES),
        expand(quast_tsv("{isolate}"), isolate=SAMPLES),
        expand(busco_txt("{isolate}"), isolate=SAMPLES),
        expand(busco_html("{isolate}"), isolate=SAMPLES),
        expand(prokka_gff("{isolate}"), isolate=SAMPLES),
        expand(prokka_faa("{isolate}"), isolate=SAMPLES),
        expand(multiqc_html("{isolate}"), isolate=SAMPLES),
        expand(multiqc_custom_tsv("{isolate}"), isolate=SAMPLES)

###############################################################################
# SECTION 5 — INPUT DISCOVERY
###############################################################################

def illumina_r1_files(i):
    if i in LIBRARY_SAMPLES:
        return [g["r1"] for lib in LIBRARY_SAMPLES[i]["libraries"] for g in lib["read_groups"]]
    d = raw_illumina_dir(i)
    pats = ["*_R1*.fastq.gz", "*_R1*.fq.gz", "*_R1*.fastq", "*_R1*.fq"]
    return sorted(sum((glob.glob(os.path.join(d, p)) for p in pats), []))

def illumina_r2_files(i):
    if i in LIBRARY_SAMPLES:
        return [g["r2"] for lib in LIBRARY_SAMPLES[i]["libraries"] for g in lib["read_groups"]]
    d = raw_illumina_dir(i)
    pats = ["*_R2*.fastq.gz", "*_R2*.fq.gz", "*_R2*.fastq", "*_R2*.fq"]
    return sorted(sum((glob.glob(os.path.join(d, p)) for p in pats), []))

def paired_illumina_files(i):
    if i in LIBRARY_SAMPLES: return illumina_r1_files(i), illumina_r2_files(i)
    return pair_illumina_files(illumina_r1_files(i), illumina_r2_files(i))

def nanopore_files(i):
    if i in LIBRARY_SAMPLES: return []
    d = raw_nanopore_dir(i)
    pats = ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]
    return sorted(sum((glob.glob(os.path.join(d, p)) for p in pats), []))

def qjoin(paths):
    return " ".join(shlex.quote(p) for p in paths)

###############################################################################
# SECTION 6 — PREFLIGHT VALIDATION AND MODE RESOLUTION
###############################################################################

RUN_MODE = {}
ISO_OVERRIDES = {}

if not SETUP_ONLY and not ADMIN_MODE:
    for iso in SAMPLES:
        r1 = illumina_r1_files(iso)
        r2 = illumina_r2_files(iso)
        r1, r2 = paired_illumina_files(iso)
        validate_input_files(r1 + r2 + nanopore_files(iso))
        ont = nanopore_files(iso)
        has_illumina = len(r1) > 0 and len(r2) > 0
        has_ont = len(ont) > 0

        if has_ont and has_illumina:
            mode = "hybrid"
        elif has_ont and not has_illumina:
            if not ALLOW_ONT_ONLY:
                raise ValueError(f"{iso}: ONT-only input detected but allow_ont_only=false")
            mode = "ont_only"
        elif has_illumina and not has_ont:
            if not ALLOW_ILLUMINA_ONLY:
                raise ValueError(f"{iso}: Illumina-only input detected but allow_illumina_only=false")
            mode = "illumina_only"
        else:
            raise ValueError(f"{iso}: no valid inputs found in raw/illumina and raw/nanopore")

        requested_mode = config["run_mode"]
        if requested_mode != "auto" and requested_mode != mode:
            raise ValueError(f"{iso}: explicit run_mode={requested_mode!r} conflicts with detected {mode!r}; no inputs will be silently omitted")

        taxid = normalize_taxid(SAMPLE_INFO[iso]["target_taxid"] or DEFAULT_TARGET_TAXID)
        lineage = SAMPLE_INFO[iso]["busco_lineage"] or DEFAULT_BUSCO_LINEAGE
        medaka_model = SAMPLE_INFO[iso]["medaka_model"] or DEFAULT_MEDAKA_MODEL

        if KRAKEN_DECISION_MODE == "manual" and REQUIRE_TAXID_FOR_KRAKEN and not taxid:
            raise ValueError(f"{iso}: Kraken manual mode requires target_taxid but none was provided")

        if RUN_BUSCO and not busco_db_ready(lineage):
            raise ValueError(
                f"{iso}: BUSCO lineage database '{lineage}' missing at {busco_lineage_dir(lineage)}. "
                f"Run: snakemake -s Hybas.smk download_busco_db --use-conda --config admin_mode=true"
            )

        RUN_MODE[iso] = mode
        ISO_OVERRIDES[iso] = {
            "target_taxid": taxid,
            "busco_lineage": lineage,
            "medaka_model": medaka_model,
        }

    budget = workflow.global_resources.get("mem_mb")
    # Newer Snakemake exposes Resource objects; older releases use scalars.
    budget = getattr(budget, "value", budget)
    reservations = {"assembly": ASSEMBLY_MEM_MB}
    if KRAKEN_PROFILE_ENABLED: reservations["kraken"] = KRAKEN_MEM_MB
    if any(planned_stages(i) for i in SAMPLES): reservations["polishing"] = POLISH_MEM_MB
    if RUN_BUSCO: reservations["busco"] = rule_memory("busco", 4096)
    if RUN_PROKKA: reservations["prokka"] = rule_memory("prokka", 4096)
    if budget is not None:
        for stage, needed in reservations.items():
            if needed > int(budget):
                raise ValueError(f"{stage} requires a {needed} MB reservation but the scheduling budget is {budget} MB. "
                    "Set total_mem_mb or --resources mem_mb to a budget the host can supply, "
                    "or reduce the stage's requirement. Snakemake must not silently cap this reservation.")
    if not os.access(BASE_DIR, os.W_OK):
        raise ValueError(f"base_dir is not writable: {BASE_DIR}")

    if KRAKEN_PROFILE_ENABLED and not kraken_db_ready():
        raise ValueError(
            f"Kraken database missing or incomplete at {KRAKEN_DB_DIR}. "
            f"Run: snakemake -s Hybas.smk download_kraken_db --use-conda --config admin_mode=true"
        )

else:
    for iso in SAMPLES:
        RUN_MODE[iso] = "setup_only" if SETUP_ONLY else "admin_only"
        ISO_OVERRIDES[iso] = {
            "target_taxid": normalize_taxid(SAMPLE_INFO[iso]["target_taxid"] or DEFAULT_TARGET_TAXID),
            "busco_lineage": SAMPLE_INFO[iso]["busco_lineage"] or DEFAULT_BUSCO_LINEAGE,
            "medaka_model": SAMPLE_INFO[iso]["medaka_model"] or DEFAULT_MEDAKA_MODEL,
        }

SETTINGS.set_sample_plans({iso: {
    "requested_mode": config["run_mode"], "mode": RUN_MODE[iso],
    "mode_reason": "Setup/admin mode" if SETUP_ONLY or ADMIN_MODE else (
        "Detected input layout" if config["run_mode"] == "auto" else "Explicit mode verified against inputs"),
    "inputs": [os.path.abspath(p) for p in illumina_r1_files(iso) + illumina_r2_files(iso) + nanopore_files(iso)],
    "overrides": ISO_OVERRIDES[iso], "polish_mode": POLISH_MODE,
    **({"library_plan": LIBRARY_SAMPLES[iso], "assembler": "spades"} if iso in LIBRARY_SAMPLES else {}),
    "planned_polishing_stages": planned_stages(iso),
    "modules": {"kraken": KRAKEN_PROFILE_ENABLED, "busco": RUN_BUSCO,
                "prokka": RUN_PROKKA, "multiqc": RUN_MULTIQC,
                "polypolish": needs_polypolish(iso)},
} for iso in SAMPLES})

###############################################################################
# SECTION 7 — STEP LABELING
###############################################################################

def logmsg(step, text):
    return text

###############################################################################
# GLOBAL PROVENANCE RULES
#
# These run once per workflow, not once per isolate.
###############################################################################

rule save_condarc:
    input:
        [PROJECT_CONDARC] if os.path.isfile(PROJECT_CONDARC) else []
    params:
        source=PROJECT_CONDARC
    output:
        provenance_condarc_file()
    message:
        "Saving project-scoped .condarc used for this HybAs run"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "save_condarc.tsv")
    run:
        os.makedirs(provenance_conda_dir(), exist_ok=True)
        if os.path.exists(PROJECT_CONDARC):
            import shutil
            shutil.copy(PROJECT_CONDARC, output[0])
        else:
            with open(output[0], "w", encoding="utf-8") as out:
                out.write("# No project .condarc found in workflow root\n")

rule setup_envs:
    input:
        lambda wc: environment_outputs(),
        provenance_condarc_file(),
        provenance_system_info(),
        provenance_workflow_info(),
        provenance_input_hashes(),


rule global_system_info:
    output:
        provenance_system_info()
    message:
        "Capturing global system/runtime provenance"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "global_system_info.tsv")
    run:
        import subprocess
        os.makedirs(provenance_runtime_dir(), exist_ok=True)

        def cmd(text):
            try:
                return subprocess.check_output(text, shell=True, text=True, stderr=subprocess.STDOUT).strip()
            except Exception:
                return "NA"

        rows = [
            ("timestamp", datetime.now().isoformat(timespec="seconds")),
            ("platform_system", platform.system()),
            ("platform_release", platform.release()),
            ("platform_version", platform.version()),
            ("platform_machine", platform.machine()),
            ("python_runtime", platform.python_version()),
            ("hostname", cmd("hostname")),
            ("uname", cmd("uname -a")),
            ("cpu_info", cmd("lscpu | sed -n '1,20p'")),
            ("memory_info", cmd("free -h")),
            ("conda_version", cmd("conda --version")),
            ("mamba_version", cmd("mamba --version")),
            ("snakemake_version", cmd("snakemake --version")),
        ]

        with open(output[0], "w", encoding="utf-8") as out:
            out.write("field\tvalue\n")
            for k, v in rows:
                out.write(f"{k}\t{str(v).replace(chr(10), ' | ')}\n")

rule global_workflow_provenance:
    params:
        effective_config=EFFECTIVE_CONFIG, version=VERSION
    input:
        save_condarc=provenance_condarc_file(), config_files=CONFIG_FILES, samples=SAMPLE_FILES, source=SOURCE_FILES
    output:
        provenance_workflow_info()
    message:
        "Capturing workflow/config/sample-sheet provenance"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "global_workflow_provenance.tsv")
    run:
        os.makedirs(provenance_runtime_dir(), exist_ok=True)

        # Resolve correct config file (works for CI and normal runs)
        try:
            cfg_path = CONFIG_FILES[-1] if CONFIG_FILES else "config.yaml"
        except (AttributeError, IndexError):
            cfg_path = "config.yaml"

        cfg_checksum = sha256_file(cfg_path) if os.path.exists(cfg_path) else "NA"
        sample_checksum = sha256_file(SAMPLE_SHEET) if SAMPLE_SHEET and os.path.exists(SAMPLE_SHEET) else "NA"
        condarc_checksum = sha256_file(input.save_condarc) if os.path.exists(input.save_condarc) else "NA"

        with open(output[0], "w", encoding="utf-8") as out:
            out.write("field\tvalue\n")
            out.write(f"workflow_version\t{WORKFLOW_VERSION}\n")
            out.write(f"source_sha256\t{sha256_file(SOURCE_FILES[0])}\n")
            out.write(f"effective_config_json\t{EFFECTIVE_CONFIG}\n")
            out.write(f"requested_config_json\t{REQUESTED_CONFIG}\n")
            out.write(f"resolved_settings_id\t{SETTINGS.settings_id}\n")
            out.write(f"resolved_settings_path\t{SETTINGS.path}\n")
            out.write(f"base_dir\t{BASE_DIR}\n")
            out.write(f"sample_sheet\t{SAMPLE_SHEET}\n")
            out.write(f"threads_default\t{THREADS}\n")
            out.write(f"keep_percent\t{KEEP_PERCENT}\n")
            out.write(f"racon_rounds\t{RACON_ROUNDS}\n")
            out.write(f"polish_mode\t{POLISH_MODE}\n")
            out.write(f"ci_mode\t{CI_MODE}\n")
            out.write(f"run_kraken\t{RUN_KRAKEN}\n")
            out.write(f"run_busco\t{RUN_BUSCO}\n")
            out.write(f"run_prokka\t{RUN_PROKKA}\n")
            out.write(f"run_multiqc\t{RUN_MULTIQC}\n")
            out.write(f"run_polypolish\t{RUN_POLYPOLISH}\n")
            out.write(f"db_root\t{DB_ROOT}\n")
            out.write(f"kraken_db_dir\t{KRAKEN_DB_DIR}\n")
            out.write(f"busco_db_dir\t{BUSCO_DB_DIR}\n")
            out.write(f"config_path\t{cfg_path}\n")
            out.write(f"config_version_label\t{CONFIG_VERSION_LABEL}\n")
            out.write(f"config_sha256\t{cfg_checksum}\n")
            out.write(f"sample_sheet_sha256\t{sample_checksum}\n")
            out.write(f"condarc_sha256\t{condarc_checksum}\n")

rule capture_env_assembly_global:
    output:
        env_conda_list_global("assembly"),
        env_export_global("assembly"),
        env_explicit_global("assembly"),
        env_pip_freeze_global("assembly"),
        env_python_version_global("assembly")
    params:
        outdir=provenance_conda_dir()
    conda: "workflow/envs/assembly.yml"
    message:
        "Capturing assembly environment manifest"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "capture_env_assembly_global.tsv")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        conda list > "{output[0]}"
        conda env export --no-builds > "{output[1]}"
        conda list --explicit > "{output[2]}"
        python -m pip freeze > "{output[3]}" || : > "{output[3]}"
        python -V > "{output[4]}" 2>&1
        """

rule capture_env_library_assembly_global:
    output:
        env_conda_list_global("library_assembly"),
        env_export_global("library_assembly"),
        env_explicit_global("library_assembly"),
        env_pip_freeze_global("library_assembly"),
        env_python_version_global("library_assembly")
    params:
        outdir=provenance_conda_dir()
    conda: "workflow/envs/library_assembly.yml"
    message:
        "Capturing library_assembly environment manifest"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "capture_env_library_assembly_global.tsv")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        conda list > "{output[0]}"
        conda env export --no-builds > "{output[1]}"
        conda list --explicit > "{output[2]}"
        python -m pip freeze > "{output[3]}" || : > "{output[3]}"
        python -V > "{output[4]}" 2>&1
        """

rule capture_env_qc_global:
    output:
        env_conda_list_global("qc"),
        env_export_global("qc"),
        env_explicit_global("qc"),
        env_pip_freeze_global("qc"),
        env_python_version_global("qc")
    params:
        outdir=provenance_conda_dir()
    conda: "workflow/envs/qc.yml"
    message:
        "Capturing QC environment manifest"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "capture_env_qc_global.tsv")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        conda list > "{output[0]}"
        conda env export --no-builds > "{output[1]}"
        conda list --explicit > "{output[2]}"
        python -m pip freeze > "{output[3]}" || : > "{output[3]}"
        python -V > "{output[4]}" 2>&1
        """

rule capture_env_polishing_global:
    output:
        env_conda_list_global("polishing"),
        env_export_global("polishing"),
        env_explicit_global("polishing"),
        env_pip_freeze_global("polishing"),
        env_python_version_global("polishing")
    params:
        outdir=provenance_conda_dir()
    conda: "workflow/envs/polishing.yml"
    message:
        "Capturing polishing environment manifest"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "capture_env_polishing_global.tsv")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        conda list > "{output[0]}"
        conda env export --no-builds > "{output[1]}"
        conda list --explicit > "{output[2]}"
        python -m pip freeze > "{output[3]}" || : > "{output[3]}"
        python -V > "{output[4]}" 2>&1
        """

rule capture_env_medaka_global:
    output:
        env_conda_list_global("medaka"),
        env_export_global("medaka"),
        env_explicit_global("medaka"),
        env_pip_freeze_global("medaka"),
        env_python_version_global("medaka")
    params:
        outdir=provenance_conda_dir()
    conda: "workflow/envs/medaka.yml"
    message:
        "Capturing medaka environment manifest"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "capture_env_medaka_global.tsv")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        conda list > "{output[0]}"
        conda env export --no-builds > "{output[1]}"
        conda list --explicit > "{output[2]}"
        python -m pip freeze > "{output[3]}" || : > "{output[3]}"
        python -V > "{output[4]}" 2>&1
        """

rule capture_env_annotation_global:
    output:
        env_conda_list_global("annotation"),
        env_export_global("annotation"),
        env_explicit_global("annotation"),
        env_pip_freeze_global("annotation"),
        env_python_version_global("annotation")
    params:
        outdir=provenance_conda_dir()
    conda: "workflow/envs/annotation.yml"
    message:
        "Capturing annotation environment manifest"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "capture_env_annotation_global.tsv")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        conda list > "{output[0]}"
        conda env export --no-builds > "{output[1]}"
        conda list --explicit > "{output[2]}"
        python -m pip freeze > "{output[3]}" || : > "{output[3]}"
        python -V > "{output[4]}" 2>&1
        """

rule database_hashes:
    input:
        kraken=lambda wc: DB_CORE if KRAKEN_PROFILE_ENABLED else [],
        busco=lambda wc: [os.path.join(busco_lineage_dir(ISO_OVERRIDES[i]["busco_lineage"]), "dataset.cfg") for i in SAMPLES] if RUN_BUSCO and not (SETUP_ONLY or ADMIN_MODE) else []
    params:
        kraken_dir=KRAKEN_DB_DIR, busco_dir=BUSCO_DB_DIR, enabled=KRAKEN_PROFILE_ENABLED
    output:
        provenance_db_hashes()
    message:
        "Hashing database core files for provenance"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "database_hashes.tsv")
    run:
        os.makedirs(provenance_db_dir(), exist_ok=True)

        rows = []
        kraken_core = [
            os.path.join(KRAKEN_DB_DIR, "hash.k2d"),
            os.path.join(KRAKEN_DB_DIR, "opts.k2d"),
            os.path.join(KRAKEN_DB_DIR, "taxo.k2d"),
        ]
        for p in kraken_core:
            rows.append(("kraken", p, sha256_file(p) if os.path.exists(p) else "NA"))

        seen_lineages = sorted(set(ISO_OVERRIDES[i]["busco_lineage"] for i in SAMPLES))
        for lineage in seen_lineages:
            lineage_path = busco_lineage_dir(lineage)
            if os.path.isdir(lineage_path):
                rows.append(("busco_lineage_dir", lineage_path, "DIR"))
                marker = os.path.join(lineage_path, ".complete")
                rows.append(("busco_lineage_complete_marker", marker, sha256_file(marker) if os.path.exists(marker) else "NA"))
            else:
                rows.append(("busco_lineage_dir", lineage_path, "NA"))

        with open(output[0], "w", encoding="utf-8") as out:
            out.write("database_type\tpath\tsha256\n")
            for dbtype, path, digest in rows:
                out.write(f"{dbtype}\t{path}\t{digest}\n")

rule global_input_hashes:
    input:
        config_files=CONFIG_FILES, samples=SAMPLE_FILES, source=SOURCE_FILES
    params:
        effective_config=EFFECTIVE_CONFIG
    output:
        provenance_input_hashes()
    message:
        "Hashing global workflow inputs"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "global_input_hashes.tsv")
    run:
        os.makedirs(provenance_runtime_dir(), exist_ok=True)

        rows = []
        for path in CONFIG_FILES:
            rows.append(("config_file", path, sha256_file(path)))
        for path in SOURCE_FILES:
            rows.append(("workflow_source", path, sha256_file(path)))
        if SAMPLE_SHEET and os.path.exists(SAMPLE_SHEET):
            rows.append(("sample_sheet", os.path.abspath(SAMPLE_SHEET), sha256_file(SAMPLE_SHEET)))
        if os.path.exists(PROJECT_CONDARC):
            rows.append((".condarc", os.path.abspath(PROJECT_CONDARC), sha256_file(PROJECT_CONDARC)))

        with open(output[0], "w", encoding="utf-8") as out:
            out.write("label\tpath\tsha256\n")
            for label, path, digest in rows:
                out.write(f"{label}\t{path}\t{digest}\n")

###############################################################################
# SECTION 8 — DATABASE ADMIN TARGETS
###############################################################################

rule download_kraken_db:
    output:
        touch(os.path.join(KRAKEN_DB_DIR, ".complete"))
    params:
        ready=kraken_db_ready(), db_identity=[file_identity(x) for x in DB_CORE],
        dbdir=KRAKEN_DB_DIR,
        hashf=os.path.join(KRAKEN_DB_DIR, "hash.k2d"),
        optsf=os.path.join(KRAKEN_DB_DIR, "opts.k2d"),
        taxof=os.path.join(KRAKEN_DB_DIR, "taxo.k2d")
    conda: "workflow/envs/qc.yml"
    log: os.path.join(BASE_DIR, "provenance", "databases", "kraken_build.log")
    resources: mem_mb=KRAKEN_MEM_MB, kraken_jobs=1
    threads: THREADS
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "download_kraken_db.tsv")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        exec > "{log}" 2>&1
        set -euo pipefail
        mkdir -p "{params.dbdir}"
        if [ -s "{params.hashf}" ] && [ -s "{params.optsf}" ] && [ -s "{params.taxof}" ]; then
          echo "Kraken DB already present at {params.dbdir}; skipping."
        else
          kraken2-build --standard --db "{params.dbdir}" --threads {threads}
        fi
        test -s "{params.hashf}" && test -s "{params.optsf}" && test -s "{params.taxof}"
        touch "{output}"
        """

rule download_busco_lineage:
    input:
        checker=os.path.join(SOURCE_ROOT, "workflow", "scripts", "check_database.py"), support=SUPPORT_FILE
    output:
        touch(os.path.join(BUSCO_DB_DIR, "{lineage}", ".complete"))
    params:
        dbdir=BUSCO_DB_DIR,
        lineage_dir=lambda wc: os.path.join(BUSCO_DB_DIR, wc.lineage),
        checker=os.path.join(SOURCE_ROOT, "workflow", "scripts", "check_database.py"),
        readiness=lambda wc: valid_busco_db(os.path.join(BUSCO_DB_DIR, wc.lineage)),
        dataset=lambda wc: file_identity(os.path.join(BUSCO_DB_DIR, wc.lineage, "dataset.cfg"))
    conda: "workflow/envs/qc.yml"
    log: os.path.join(BASE_DIR, "provenance", "databases", "busco_{lineage}.log")
    threads: 1
    message: "Prepare and verify BUSCO lineage {wildcards.lineage}"
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "download_busco_lineage_{lineage}.tsv")
    shell:
        r"""
        mkdir -p "{params.dbdir}" "$(dirname "{log}")"
        (
          if python "{params.checker}" "{params.lineage_dir}"; then
            echo "Verified BUSCO dataset already present"
          else
            tempdir=$(mktemp -d "{params.dbdir}/.download.{wildcards.lineage}.XXXXXX")
            busco --download "{wildcards.lineage}" --download_path "$tempdir"
            python "{params.checker}" "$tempdir/{wildcards.lineage}"
            if [ -e "{params.lineage_dir}" ]; then
              mv "{params.lineage_dir}" "$tempdir/previous-incomplete"
            fi
            mv "$tempdir/{wildcards.lineage}" "{params.lineage_dir}"
          fi
          touch "{output}"
        ) > "{log}" 2>&1
        """

rule download_busco_db:
    input:
        expand(os.path.join(BUSCO_DB_DIR, "{lineage}", ".complete"), lineage=BUSCO_DOWNLOAD_LINEAGES)

rule update_databases:
    input:
        os.path.join(KRAKEN_DB_DIR, ".complete"),
        expand(os.path.join(BUSCO_DB_DIR, "{lineage}", ".complete"), lineage=BUSCO_DOWNLOAD_LINEAGES)


###############################################################################
# SECTION 10 — SYSTEM AND WORKFLOW PROVENANCE (PER ISOLATE)
###############################################################################

rule system_info:
    output:
        report(
            system_info_tsv("{isolate}"),
            caption="workflow/report/system_info_tsv.rst",
            category="Run metadata",
            subcategory="System info",
        )
    benchmark:
        bench("{isolate}", "system_info")
    run:
        os.makedirs(metadata_dir(wildcards.isolate), exist_ok=True)
        with open(output[0], "w", encoding="utf-8") as out:
            out.write("field\tvalue\n")
            out.write(f"timestamp\t{datetime.now().isoformat(timespec='seconds')}\n")
            out.write(f"platform_system\t{platform.system()}\n")
            out.write(f"platform_release\t{platform.release()}\n")
            out.write(f"platform_version\t{platform.version()}\n")
            out.write(f"platform_machine\t{platform.machine()}\n")
            out.write(f"python_runtime\t{platform.python_version()}\n")

rule workflow_provenance:
    input:
        global_info=provenance_workflow_info(), config_files=CONFIG_FILES, samples=SAMPLE_FILES
    params:
        effective_config=EFFECTIVE_CONFIG
    output:
        report(
            workflow_provenance_tsv("{isolate}"),
            caption="workflow/report/workflow_provenance_tsv.rst",
            category="Run metadata",
            subcategory="Workflow provenance",
        )
    benchmark:
        bench("{isolate}", "workflow_provenance")
    run:
        os.makedirs(metadata_dir(wildcards.isolate), exist_ok=True)
        cfg_checksum = sha256_file(CONFIG_FILES[-1]) if CONFIG_FILES else "NA"
        sample_checksum = sha256_file(SAMPLE_SHEET) if SAMPLE_SHEET and os.path.exists(SAMPLE_SHEET) else "NA"

        with open(output[0], "w", encoding="utf-8") as out:
            out.write("field\tvalue\n")
            out.write(f"workflow_version\t{WORKFLOW_VERSION}\n")
            out.write(f"source_sha256\t{sha256_file(SOURCE_FILES[0])}\n")
            out.write(f"effective_config_json\t{EFFECTIVE_CONFIG}\n")
            out.write(f"requested_config_json\t{REQUESTED_CONFIG}\n")
            out.write(f"resolved_settings_id\t{SETTINGS.settings_id}\n")
            out.write(f"resolved_settings_path\t{SETTINGS.path}\n")
            out.write(f"threads_default\t{THREADS}\n")
            out.write(f"keep_percent\t{KEEP_PERCENT}\n")
            out.write(f"racon_rounds\t{RACON_ROUNDS}\n")
            out.write(f"polish_mode\t{POLISH_MODE}\n")
            out.write(f"run_kraken\t{RUN_KRAKEN}\n")
            out.write(f"run_busco\t{RUN_BUSCO}\n")
            out.write(f"run_prokka\t{RUN_PROKKA}\n")
            out.write(f"run_multiqc\t{RUN_MULTIQC}\n")
            out.write(f"run_polypolish\t{RUN_POLYPOLISH}\n")
            out.write(f"db_root\t{DB_ROOT}\n")
            out.write(f"kraken_db_dir\t{KRAKEN_DB_DIR}\n")
            out.write(f"busco_db_dir\t{BUSCO_DB_DIR}\n")
            out.write(f"config_sha256\t{cfg_checksum}\n")
            out.write(f"sample_sheet_sha256\t{sample_checksum}\n")

###############################################################################
# SECTION 11 — RAW MERGING
###############################################################################


rule merge_illumina_libraries:
    wildcard_constraints:
        isolate=LIBRARY_PATTERN
    input:
        r1=lambda wc: illumina_r1_files(wc.isolate),
        r2=lambda wc: illumina_r2_files(wc.isolate),
        adapter=LIBRARY_ADAPTER
    output:
        r1=merged_r1("{isolate}"), r2=merged_r2("{isolate}"),
        audit=directory(library_audit("{isolate}"))
    params:
        plan=lambda wc: json.dumps(LIBRARY_SAMPLES[wc.isolate], sort_keys=True),
        outdir=lambda wc: work_inputs_merged_dir(wc.isolate)
    log: os.path.join(logs_dir("{isolate}"), "merge_illumina_libraries.log")
    benchmark: bench("{isolate}", "merge_illumina_libraries")
    conda: "workflow/envs/assembly.yml"
    threads: 1
    message: "Audit read pairs and preserve library identities: {wildcards.isolate}"
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        python {input.adapter:q} prepare --plan {params.plan:q} --r1 {output.r1:q} --r2 {output.r2:q} --audit {output.audit:q} > {log:q} 2>&1
        """

rule fastp_libraries:
    wildcard_constraints:
        isolate=LIBRARY_PATTERN
    input:
        r1=clean_r1("{isolate}"), r2=clean_r2("{isolate}"),
        audit=library_audit("{isolate}"), adapter=LIBRARY_ADAPTER
    output:
        r1=trim_r1("{isolate}"), r2=trim_r2("{isolate}"),
        html=report(os.path.join(fastp_dir("{isolate}"), "fastp.html"), caption="workflow/report/library_fastp.rst", category="QC", subcategory="Library trimming"),
        json=report(os.path.join(fastp_dir("{isolate}"), "fastp.json"), caption="workflow/report/library_fastp.rst", category="QC", subcategory="Library trimming"),
        libraries=directory(library_trimmed("{isolate}"))
    params:
        trimdir=lambda wc: trimmed_dir(wc.isolate),
        reportdir=lambda wc: fastp_dir(wc.isolate)
    log: os.path.join(logs_dir("{isolate}"), "fastp_libraries.log")
    benchmark: bench("{isolate}", "fastp_libraries")
    conda: "workflow/envs/assembly.yml"
    threads: THREADS
    message: "Trim each Illumina library separately: {wildcards.isolate}"
    shell:
        r"""
        mkdir -p "{params.trimdir}" "{params.reportdir}" "$(dirname "{log}")"
        python {input.adapter:q} trim --plan "{input.audit}/plan.json" --r1 {input.r1:q} --r2 {input.r2:q} \
          --outdir {output.libraries:q} --combined1 {output.r1:q} --combined2 {output.r2:q} \
          --report-html {output.html:q} --report-json {output.json:q} --threads {threads} > {log:q} 2>&1
        """

rule spades_libraries:
    wildcard_constraints:
        isolate=LIBRARY_PATTERN
    input:
        libraries=library_trimmed("{isolate}"), adapter=LIBRARY_ADAPTER
    output:
        bundle=directory(spades_outdir("{isolate}")), snapshot=stage_spades("{isolate}")
    params:
        parent=lambda wc: work_assembly_dir(wc.isolate),
        memory_gb=SPADES_MEMORY_GB, minimum=MIN_FASTA_LENGTH
    log: os.path.join(logs_dir("{isolate}"), "spades_libraries.log")
    benchmark: bench("{isolate}", "spades_libraries")
    conda: "workflow/envs/library_assembly.yml"
    resources:
        mem_mb=ASSEMBLY_MEM_MB
    threads: THREADS
    message: "SPAdes assembly with separate paired-end libraries: {wildcards.isolate}"
    shell:
        r"""
        mkdir -p "{params.parent}" "$(dirname "{log}")"
        (
          set -euo pipefail
          tempdir=$(mktemp -d "{params.parent}/.spades.XXXXXX")
          echo "SPAdes work directory: $tempdir"
          spades.py --dataset "{input.libraries}/dataset.yaml" -o "$tempdir" -t {threads} -m {params.memory_gb} --phred-offset 33
          mv "$tempdir" "{output.bundle}"
          python {input.adapter:q} promote "{output.bundle}/contigs.fasta" {output.snapshot:q} {params.minimum}
        ) > {log:q} 2>&1
        """

rule merge_illumina:
    wildcard_constraints:
        isolate=LEGACY_PATTERN
    input:
        r1=lambda wc: paired_illumina_files(wc.isolate)[0],
        r2=lambda wc: paired_illumina_files(wc.isolate)[1]
    output:
        r1=merged_r1("{isolate}"),
        r2=merged_r2("{isolate}")
    params:
        r1_list=lambda wc, input: qjoin(input.r1),
        r2_list=lambda wc, input: qjoin(input.r2),
        outdir=lambda wc: work_inputs_merged_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "merge_illumina.log")
    benchmark:
        bench("{isolate}", "merge_illumina")
    conda:
        "workflow/envs/assembly.yml"
    threads: 1
    message:
        logmsg(1, "Merging Illumina lanes")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"
        if [ -z "{params.r1_list}" ] || [ -z "{params.r2_list}" ]; then
          : > "{output.r1}"
          : > "{output.r2}"
          exit 0
        fi
        tmp1="{params.outdir}/.r1.tmp.fq"
        tmp2="{params.outdir}/.r2.tmp.fq"
        > "$tmp1"; > "$tmp2"
        for f in {params.r1_list}; do case "$f" in *.gz) zcat "$f" ;; *) cat "$f" ;; esac >> "$tmp1"; done
        for f in {params.r2_list}; do case "$f" in *.gz) zcat "$f" ;; *) cat "$f" ;; esac >> "$tmp2"; done
        pigz -p {threads} -n -c "$tmp1" > "{output.r1}"
        pigz -p {threads} -n -c "$tmp2" > "{output.r2}"
        rm -f "$tmp1" "$tmp2"
        ) > "{log}" 2>&1
        """

rule merge_ont:
    input:
        files=lambda wc: nanopore_files(wc.isolate)
    output:
        merged_ont("{isolate}")
    params:
        ont_list=lambda wc, input: qjoin(input.files),
        outdir=lambda wc: work_inputs_merged_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "merge_ont.log")
    benchmark:
        bench("{isolate}", "merge_ont")
    conda:
        "workflow/envs/assembly.yml"
    threads: 1
    message:
        logmsg(2, "Merging ONT FASTQ")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"
        if [ -z "{params.ont_list}" ]; then
          : > "{output}"
          exit 0
        fi
        tmp="{params.outdir}/.ont.tmp.fastq"
        > "$tmp"
        for f in {params.ont_list}; do case "$f" in *.gz) zcat "$f" ;; *) cat "$f" ;; esac >> "$tmp"; done
        pigz -p {threads} -n -c "$tmp" > "{output}"
        rm -f "$tmp"
        ) > "{log}" 2>&1
        """

###############################################################################
# SECTION 12 — RAW QC
###############################################################################

rule fastqc_pre:
    input:
        r1=merged_r1("{isolate}"),
        r2=merged_r2("{isolate}")
    output:
        directory(fastqc_pre_dir("{isolate}"))
    params:
        outdir=lambda wc: fastqc_pre_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "fastqc_pre.log")
    conda:
        "workflow/envs/qc.yml"
    threads: 2
    message:
        logmsg(3, "FastQC raw merged Illumina")
    benchmark:
        bench("{isolate}", "fastqc_pre")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"
        if [ -s "{input.r1}" ]; then fastqc -t {threads} -o "{params.outdir}" "{input.r1}"; fi
        if [ -s "{input.r2}" ]; then fastqc -t {threads} -o "{params.outdir}" "{input.r2}"; fi
        ) > "{log}" 2>&1
        """

rule nanoplot_raw:
    input:
        merged_ont("{isolate}")
    output:
        directory(nanoplot_raw_dir("{isolate}"))
    params:
        outdir=lambda wc: nanoplot_raw_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "nanoplot_raw.log")
    conda:
        "workflow/envs/qc.yml"
    threads: 2
    message:
        logmsg(4, "NanoPlot raw ONT")
    benchmark:
        bench("{isolate}", "nanoplot_raw")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"
        if [ -s "{input}" ]; then
          NanoPlot \
            --fastq "{input}" \
            -o "{params.outdir}" \
            --threads {threads} \
            --tsv_stats \
            --prefix raw_ont_
        else
          echo "No ONT input" > "{params.outdir}/README.txt"
        fi
        ) > "{log}" 2>&1
        """

rule kraken2_illumina:
    input:
        r1=merged_r1("{isolate}"),
        r2=merged_r2("{isolate}")
    output:
        kraken=report(
            os.path.join(kraken_dir("{isolate}"), "illumina.kraken2.txt"),
            caption="workflow/report/kraken_report.rst",
            category="QC",
            subcategory="Kraken2",
        ),
        report_txt=report(
            os.path.join(kraken_dir("{isolate}"), "illumina.report.txt"),
            caption="workflow/report/kraken_report.rst",
            category="QC",
            subcategory="Kraken2",
        )
    params:
        enabled=KRAKEN_PROFILE_ENABLED, database=KRAKEN_DB_DIR,
        db_identity=lambda wc: [file_identity(x) for x in DB_CORE] if KRAKEN_PROFILE_ENABLED else [],
        outdir=lambda wc: kraken_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "kraken2_illumina.log")
    benchmark:
        bench("{isolate}", "kraken2_illumina")
    conda:
        "workflow/envs/qc.yml"
    resources:
        mem_mb=KRAKEN_MEM_MB, kraken_jobs=1
    threads: THREADS
    message:
        logmsg(5, "Kraken2 raw Illumina")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"
        if [ "{KRAKEN_PROFILE_ENABLED}" != "True" ]; then
          : > "{output.kraken}"
          : > "{output.report_txt}"
          exit 0
        fi
        if [ -s "{input.r1}" ] && [ -s "{input.r2}" ]; then
          kraken2 --db "{KRAKEN_DB_DIR}" --threads {threads} \
            --paired \
            --report "{output.report_txt}" \
            "{input.r1}" "{input.r2}" > "{output.kraken}"
        else
          : > "{output.kraken}"
          : > "{output.report_txt}"
        fi
        ) > "{log}" 2>&1
        """

rule kraken2_nanopore:
    input:
        ont=merged_ont("{isolate}")
    output:
        kraken=report(
            os.path.join(kraken_dir("{isolate}"), "nanopore.kraken2.txt"),
            caption="workflow/report/kraken_report.rst",
            category="QC",
            subcategory="Kraken2",
        ),
        report_txt=report(
            os.path.join(kraken_dir("{isolate}"), "nanopore.report.txt"),
            caption="workflow/report/kraken_report.rst",
            category="QC",
            subcategory="Kraken2",
        )
    params:
        enabled=KRAKEN_PROFILE_ENABLED, database=KRAKEN_DB_DIR,
        db_identity=lambda wc: [file_identity(x) for x in DB_CORE] if KRAKEN_PROFILE_ENABLED else [],
        outdir=lambda wc: kraken_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "kraken2_nanopore.log")
    benchmark:
        bench("{isolate}", "kraken2_nanopore")
    conda:
        "workflow/envs/qc.yml"
    resources:
        mem_mb=KRAKEN_MEM_MB, kraken_jobs=1
    threads: THREADS
    message:
        logmsg(6, "Kraken2 raw ONT")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"
        if [ "{KRAKEN_PROFILE_ENABLED}" != "True" ]; then
          : > "{output.kraken}"
          : > "{output.report_txt}"
          exit 0
        fi
        if [ -s "{input.ont}" ]; then
          kraken2 --db "{KRAKEN_DB_DIR}" --threads {threads} \
            --report "{output.report_txt}" "{input.ont}" > "{output.kraken}"
        else
          : > "{output.kraken}"
          : > "{output.report_txt}"
        fi
        ) > "{log}" 2>&1
        """

rule kraken_decision:
    params:
        settings=AWARE_SETTINGS,
        controller=AWARE_CODE_HASH,
        target=lambda wc: ISO_OVERRIDES[wc.isolate]["target_taxid"]
    input:
        illumina_report=os.path.join(kraken_dir("{isolate}"), "illumina.report.txt"),
        ont_report=os.path.join(kraken_dir("{isolate}"), "nanopore.report.txt")
    output:
        report(
            kraken_decision_tsv("{isolate}"),
            caption="workflow/report/run_summary_tsv.rst",
            category="Run metadata",
            subcategory="Kraken decision",
        )
    benchmark:
        bench("{isolate}", "kraken_decision")
    run:
        iso = wildcards.isolate
        os.makedirs(metadata_dir(iso), exist_ok=True)

        manual_target_taxid = normalize_taxid(ISO_OVERRIDES[iso]["target_taxid"])

        target_info = choose_selected_kraken_target(
            mode=KRAKEN_DECISION_MODE,
            manual_target_taxid=manual_target_taxid,
            illumina_report=input.illumina_report,
            ont_report=input.ont_report,
        )

        # keep the global combined decision for reporting only
        combined_decision, combined_reason = choose_kraken_decision(
            mode=KRAKEN_DECISION_MODE,
            selected_taxid=target_info["selected_taxid"],
            selected_pct_for_decision=target_info["selected_pct_for_decision"],
        )

        # per-platform cleaning decisions
        illumina_decision, illumina_reason = choose_platform_cleaning_decision(
            mode=KRAKEN_DECISION_MODE,
            selected_taxid=target_info["selected_taxid"],
            platform_pct=target_info["illumina_selected_pct"],
            platform_name="illumina",
        )

        ont_decision, ont_reason = choose_platform_cleaning_decision(
            mode=KRAKEN_DECISION_MODE,
            selected_taxid=target_info["selected_taxid"],
            platform_pct=target_info["ont_selected_pct"],
            platform_name="ont",
        )

        confidence_info = compute_taxonomic_confidence(
            illumina_report=input.illumina_report,
            ont_report=input.ont_report,
            selected_taxid=target_info["selected_taxid"],
        )

        with open(output[0], "w", encoding="utf-8") as out:
            out.write("field\tvalue\n")
            out.write(f"isolate\t{iso}\n")
            out.write(f"kraken_decision_mode\t{KRAKEN_DECISION_MODE}\n")

            out.write(f"manual_target_taxid\t{manual_target_taxid}\n")

            out.write(f"selected_target_taxid\t{target_info['selected_taxid']}\n")
            out.write(f"selected_target_name\t{target_info['selected_name']}\n")
            out.write(f"selected_target_rank\t{target_info['selected_rank']}\n")
            out.write(f"selected_target_source\t{target_info['selected_source']}\n")
            out.write(f"selected_target_percent_for_decision\t{target_info['selected_pct_for_decision']:.4f}\n")
            out.write(f"illumina_selected_target_percent\t{target_info['illumina_selected_pct']:.4f}\n")
            out.write(f"ont_selected_target_percent\t{target_info['ont_selected_pct']:.4f}\n")

            out.write(f"inferred_target_taxid\t{target_info['inferred_taxid']}\n")
            out.write(f"inferred_target_name\t{target_info['inferred_name']}\n")
            out.write(f"inferred_target_rank\t{target_info['inferred_rank']}\n")
            out.write(f"inferred_target_percent\t{target_info['inferred_pct']:.4f}\n")
            out.write(f"target_selection_reason\t{target_info['selection_reason']}\n")

            out.write(f"runner_up_taxid\t{confidence_info['runner_up_taxid']}\n")
            out.write(f"runner_up_name\t{confidence_info['runner_up_name']}\n")
            out.write(f"runner_up_percent\t{confidence_info['runner_up_pct']:.4f}\n")
            out.write(f"selected_vs_runner_up_delta_percent\t{confidence_info['delta_pct']:.4f}\n")
            out.write(f"platform_agreement_score\t{confidence_info['platform_agreement_score']:.4f}\n")
            out.write(f"platform_agreement_reason\t{confidence_info['platform_agreement_reason']}\n")
            out.write(f"decision_confidence\t{confidence_info['decision_confidence']}\n")
            out.write(f"awareness_state\t{confidence_info['awareness_state']}\n")

            out.write(f"kraken_decision\t{combined_decision}\n")
            out.write(f"kraken_decision_reason\t{combined_reason}\n")

            out.write(f"illumina_kraken_decision\t{illumina_decision}\n")
            out.write(f"illumina_kraken_decision_reason\t{illumina_reason}\n")
            out.write(f"ont_kraken_decision\t{ont_decision}\n")
            out.write(f"ont_kraken_decision_reason\t{ont_reason}\n")

            out.write(f"kraken_include_children\t{KRAKEN_INCLUDE_CHILDREN}\n")


###############################################################################
# SECTION 13 — CLEANING
###############################################################################

rule clean_illumina:
    input:
        r1=merged_r1("{isolate}"),
        r2=merged_r2("{isolate}"),
        kr=os.path.join(kraken_dir("{isolate}"), "illumina.kraken2.txt"),
        report=os.path.join(kraken_dir("{isolate}"), "illumina.report.txt"),
        decision=kraken_decision_tsv("{isolate}")
    output:
        r1=clean_r1("{isolate}"),
        r2=clean_r2("{isolate}")
    params:
        include_children=KRAKEN_INCLUDE_CHILDREN,
        outdir=lambda wc: work_inputs_cleaned_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "clean_illumina.log")
    benchmark:
        bench("{isolate}", "clean_illumina")
    conda:
        "workflow/envs/qc.yml"
    threads: 4
    message:
        logmsg(7, "Illumina cleaning / passthrough")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"

        if [ ! -s "{input.r1}" ] || [ ! -s "{input.r2}" ]; then
          : > "{output.r1}"
          : > "{output.r2}"
          exit 0
        fi

        decision=$(awk -F '\t' '$1=="illumina_kraken_decision"{{print $2}}' "{input.decision}" | tail -n1)
        selected_taxid=$(awk -F '\t' '$1=="selected_target_taxid"{{print $2}}' "{input.decision}" | tail -n1)

        if [ "$decision" = "passthrough" ] || [ -z "$selected_taxid" ]; then
          cp -a "{input.r1}" "{output.r1}"
          cp -a "{input.r2}" "{output.r2}"
          exit 0
        fi

        tmp1="{params.outdir}/.illumina_R1.clean.fq"
        tmp2="{params.outdir}/.illumina_R2.clean.fq"

        extra_children=""
        if [ "{KRAKEN_INCLUDE_CHILDREN}" = "True" ]; then
          extra_children="--include-children"
        fi

        extract_kraken_reads.py -k "{input.kr}" -r "{input.report}" \
          -s1 "{input.r1}" -s2 "{input.r2}" \
          --taxid "$selected_taxid" $extra_children \
          -o "$tmp1" -o2 "$tmp2" --fastq-output

        pigz -p {threads} -n -f "$tmp1"
        pigz -p {threads} -n -f "$tmp2"
        mv "$tmp1.gz" "{output.r1}"
        mv "$tmp2.gz" "{output.r2}"
        gzip -t "{output.r1}"
        gzip -t "{output.r2}"
        ) > "{log}" 2>&1
        """

rule clean_ont:
    input:
        ont=merged_ont("{isolate}"),
        kr=os.path.join(kraken_dir("{isolate}"), "nanopore.kraken2.txt"),
        report=os.path.join(kraken_dir("{isolate}"), "nanopore.report.txt"),
        decision=kraken_decision_tsv("{isolate}")
    output:
        cleaned=clean_ont_path("{isolate}")
    params:
        include_children=KRAKEN_INCLUDE_CHILDREN,
        outdir=lambda wc: work_inputs_cleaned_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "clean_ont.log")
    benchmark:
        bench("{isolate}", "clean_ont")
    conda:
        "workflow/envs/qc.yml"
    threads: 4
    message:
        logmsg(8, "ONT cleaning / passthrough")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"

        if [ ! -s "{input.ont}" ]; then
          : > "{output.cleaned}"
          exit 0
        fi

        decision=$(awk -F '\t' '$1=="ont_kraken_decision"{{print $2}}' "{input.decision}" | tail -n1)
        selected_taxid=$(awk -F '\t' '$1=="selected_target_taxid"{{print $2}}' "{input.decision}" | tail -n1)

        if [ "$decision" = "passthrough" ] || [ -z "$selected_taxid" ]; then
          cp -a "{input.ont}" "{output.cleaned}"
          exit 0
        fi

        tmp="{params.outdir}/.ont.clean.fastq"

        extra_children=""
        if [ "{KRAKEN_INCLUDE_CHILDREN}" = "True" ]; then
          extra_children="--include-children"
        fi

        extract_kraken_reads.py -k "{input.kr}" -r "{input.report}" \
          -s "{input.ont}" --taxid "$selected_taxid" $extra_children \
          -o "$tmp" --fastq-output

        pigz -p {threads} -n -f "$tmp"
        mv "$tmp.gz" "{output.cleaned}"
        gzip -t "{output.cleaned}"
        ) > "{log}" 2>&1
        """

###############################################################################
# SECTION 14 — PREPROCESSING
###############################################################################

rule fastp_trim:
    wildcard_constraints:
        isolate=LEGACY_PATTERN
    input:
        r1=clean_r1("{isolate}"),
        r2=clean_r2("{isolate}")
    output:
        r1=trim_r1("{isolate}"),
        r2=trim_r2("{isolate}"),
        html=report(
            os.path.join(fastp_dir("{isolate}"), "fastp.html"),
            caption="workflow/report/fastp_html.rst",
            category="QC",
            subcategory="fastp",
        ),
        json=report(
            os.path.join(fastp_dir("{isolate}"), "fastp.json"),
            caption="workflow/report/fastp_json.rst",
            category="QC",
            subcategory="fastp",
        )
    params:
        trimdir=lambda wc: trimmed_dir(wc.isolate),
        fastpdir=lambda wc: fastp_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "fastp_trim.log")
    benchmark:
        bench("{isolate}", "fastp_trim")
    conda:
        "workflow/envs/assembly.yml"
    threads: THREADS
    message:
        logmsg(9, "fastp trimming")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.trimdir}" "{params.fastpdir}" "{params.logdir}"
        if [ -s "{input.r1}" ] && [ -s "{input.r2}" ]; then
          fastp -i "{input.r1}" -I "{input.r2}" -o "{output.r1}" -O "{output.r2}" \
            --thread {threads} --html "{output.html}" --json "{output.json}"
        else
          : > "{output.r1}"
          : > "{output.r2}"
          : > "{output.html}"
          : > "{output.json}"
        fi
        ) > "{log}" 2>&1
        """

rule filtlong_ont:
    input:
        clean_ont_path("{isolate}")
    output:
        filt_ont("{isolate}")
    params:
        keep_percent=KEEP_PERCENT,
        outdir=lambda wc: work_inputs_filtered_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "filtlong_ont.log")
    benchmark:
        bench("{isolate}", "filtlong_ont")
    conda:
        "workflow/envs/assembly.yml"
    threads: THREADS
    message:
        logmsg(10, "Filtlong ONT")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"
        if [ -s "{input}" ]; then
          filtlong --min_length 1000 --keep_percent {KEEP_PERCENT} "{input}" > "{output}.tmp.fastq"
          pigz -p {threads} -n -c "{output}.tmp.fastq" > "{output}"
          rm -f "{output}.tmp.fastq"
        else
          : > "{output}"
        fi
        ) > "{log}" 2>&1
        """

rule nanoplot_filt:
    input:
        filt_ont("{isolate}")
    output:
        directory(nanoplot_filt_dir("{isolate}"))
    params:
        outdir=lambda wc: nanoplot_filt_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "nanoplot_filt.log")
    conda:
        "workflow/envs/qc.yml"
    threads: 2
    message:
        logmsg(11, "NanoPlot filtered ONT")
    benchmark:
        bench("{isolate}", "nanoplot_filt")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"
        if [ -s "{input}" ]; then
          NanoPlot \
            --fastq "{input}" \
            -o "{params.outdir}" \
            --threads {threads} \
            --tsv_stats \
            --prefix filt_ont_
        else
          echo "No ONT input" > "{params.outdir}/README.txt"
        fi
        ) > "{log}" 2>&1
        """

###############################################################################
# SECTION 15 — MODE DECLARATION / SKIP NOTES
###############################################################################

rule run_mode:
    params:
        mode=lambda wc: RUN_MODE[wc.isolate],
        polish=POLISH_MODE,
        ci=CI_MODE,
        kraken=KRAKEN_DECISION_MODE
    output:
        run_mode_file("{isolate}")
    benchmark:
        bench("{isolate}", "run_mode")
    run:
        os.makedirs(metadata_dir(wildcards.isolate), exist_ok=True)
        with open(output[0], "w", encoding="utf-8") as out:
            out.write("field\tvalue\n")
            out.write(f"isolate\t{wildcards.isolate}\n")
            out.write(f"run_mode\t{RUN_MODE[wildcards.isolate]}\n")
            out.write(f"polish_mode_requested\t{POLISH_MODE}\n")
            out.write(f"ci_mode\t{CI_MODE}\n")
            out.write(f"kraken_decision_mode\t{KRAKEN_DECISION_MODE}\n")

rule skip_reasons:
    input:
        kd=kraken_decision_tsv("{isolate}"), execution=polish_execution("{isolate}")
    output:
        report(
            skip_reasons_tsv("{isolate}"),
            caption="workflow/report/skip_reasons_tsv.rst",
            category="Run metadata",
            subcategory="Skip reasons",
        )
    benchmark:
        bench("{isolate}", "skip_reasons")
    run:
        iso = wildcards.isolate
        mode = RUN_MODE[iso]
        rows = []
        execution = json.loads(Path(input.execution).read_text())
        rows.extend((stage, reason) for stage, reason in execution["skipped_stages"].items())

        if mode == "setup_only":
            rows.append(("workflow", "setup_only_mode"))
        if mode == "admin_only":
            rows.append(("workflow", "admin_only_mode"))

        kd = read_key_value_tsv(input.kd, key_col=0, val_col=1, has_header=True)
        decision = kd.get("kraken_decision", "NA")
        reason_detail = kd.get("kraken_decision_reason", "NA")
        selected_source = kd.get("selected_target_source", "NA")
        selected_taxid = kd.get("selected_target_taxid", "NA")
        manual_target = kd.get("manual_target_taxid", "NA")
        inferred_taxid = kd.get("inferred_target_taxid", "NA")

        awareness_state = kd.get("awareness_state", "NA")
        decision_confidence = kd.get("decision_confidence", "NA")
        runner_up_taxid = kd.get("runner_up_taxid", "NA")
        runner_up_name = kd.get("runner_up_name", "NA")
        runner_up_pct = kd.get("runner_up_percent", "NA")
        delta_pct = kd.get("selected_vs_runner_up_delta_percent", "NA")
        agreement_score = kd.get("platform_agreement_score", "NA")
        agreement_reason = kd.get("platform_agreement_reason", "NA")

        rows.append(("kraken_decision", f"{decision}; {reason_detail}"))
        rows.append(("kraken_target_selection", f"source={selected_source}; selected_taxid={selected_taxid}; manual_target={manual_target}; inferred_target={inferred_taxid}"))
        rows.append(("kraken_awareness", f"state={awareness_state}; confidence={decision_confidence}; delta_pct={delta_pct}; agreement_score={agreement_score}; agreement_reason={agreement_reason}; runner_up_taxid={runner_up_taxid}; runner_up_name={runner_up_name}; runner_up_percent={runner_up_pct}"))

        if KRAKEN_DECISION_MODE == "off":
            rows.append(("kraken", "decision_mode_off_kraken_disabled"))

        elif KRAKEN_DECISION_MODE == "manual":
            if not normalize_taxid(manual_target):
                rows.append(("kraken_whitelist", "manual_mode_empty_target_taxid_passthrough_used"))
            else:
                rows.append(("kraken_target_selection", "manual_mode_explicit_target_used"))

        elif KRAKEN_DECISION_MODE == "auto":
            if selected_source == "manual_override_in_auto":
                rows.append(("kraken_target_selection", "auto_mode_manual_override_used"))
            elif selected_source == "inferred_from_kraken":
                rows.append(("kraken_target_selection", "auto_mode_inferred_target_used"))
            elif selected_source == "auto_inference_failed":
                rows.append(("kraken_target_selection", "auto_mode_inference_failed"))

        elif KRAKEN_DECISION_MODE == "aware":
            if selected_source == "manual_override_in_aware":
                rows.append(("kraken_target_selection", "aware_mode_manual_override_used"))
            elif selected_source == "inferred_from_kraken":
                rows.append(("kraken_target_selection", "aware_mode_inferred_target_used"))
            elif selected_source == "aware_inference_failed":
                rows.append(("kraken_target_selection", "aware_mode_inference_failed"))

            rows.append(("kraken_awareness_state", awareness_state))
            rows.append(("kraken_awareness_confidence", decision_confidence))

            if awareness_state == "clean_stable":
                rows.append(("kraken_awareness_interpretation", "high-confidence taxonomically coherent sample"))
            elif awareness_state == "clean_filterable":
                rows.append(("kraken_awareness_interpretation", "moderate-confidence dominant target with recoverable contamination burden"))
            elif awareness_state == "taxonomically_ambiguous":
                rows.append(("kraken_awareness_interpretation", "low-confidence or competing taxonomic structure detected"))
            elif awareness_state == "fail_no_target":
                rows.append(("kraken_awareness_interpretation", "no meaningful target available for awareness decisioning"))

        with open(output[0], "w", encoding="utf-8") as out:
            out.write("field\tvalue\n")
            for field, value in rows:
                out.write(f"{field}\t{value}\n")

###############################################################################
# SECTION 16 — ASSEMBLY
###############################################################################

rule unicycler:
    wildcard_constraints:
        isolate=LEGACY_PATTERN
    input:
        ont=filt_ont("{isolate}"),
        r1=trim_r1("{isolate}"),
        r2=trim_r2("{isolate}"),
        mode=run_mode_file("{isolate}")
    output:
        asm=unicycler_asm("{isolate}"),
        snapshot=stage_unicycler("{isolate}")
    params:
        spades_memory_gb=SPADES_MEMORY_GB, min_fasta_length=MIN_FASTA_LENGTH,
        outdir=lambda wc: unicycler_outdir(wc.isolate),
        assemblydir=lambda wc: work_assembly_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate),
        mode=lambda wc: RUN_MODE[wc.isolate]
    log:
        os.path.join(logs_dir("{isolate}"), "unicycler.log")
    benchmark:
        bench("{isolate}", "unicycler")
    conda:
        "workflow/envs/assembly.yml"
    resources:
        mem_mb=ASSEMBLY_MEM_MB, disk_mb=rule_memory("assembly_disk", 20000)
    threads: THREADS
    message:
        logmsg(12, "Unicycler assembly")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.assemblydir}" "{params.logdir}"
        trap 'rc=$?; if [ "$rc" -ne 0 ] && [ -d "{params.outdir}" ]; then find "{params.outdir}" -type f -name "*.log" -exec cp --parents {{}} "{params.logdir}/" \; ; fi; exit "$rc"' EXIT
        if [ "{params.mode}" = "hybrid" ]; then
          unicycler -l "{input.ont}" -1 "{input.r1}" -2 "{input.r2}" -o "{params.outdir}" -t {threads} --min_fasta_length {params.min_fasta_length} --spades_options "-m {params.spades_memory_gb}"
        elif [ "{params.mode}" = "ont_only" ]; then
          unicycler -l "{input.ont}" -o "{params.outdir}" -t {threads} --min_fasta_length {params.min_fasta_length} --spades_options "-m {params.spades_memory_gb}"
        elif [ "{params.mode}" = "illumina_only" ]; then
          unicycler -1 "{input.r1}" -2 "{input.r2}" -o "{params.outdir}" -t {threads} --min_fasta_length {params.min_fasta_length} --spades_options "-m {params.spades_memory_gb}"
        else
          echo "Unknown run mode"
          exit 1
        fi
        cp "{output.asm}" "{output.snapshot}"
        ) > "{log}" 2>&1
        """

###############################################################################
# SECTION 17 — STAGE-RESOLVED POLISHING
###############################################################################

rule racon_seed:
    input:
        stage_unicycler("{isolate}")
    output:
        stage_racon0("{isolate}")
    params:
        outdir=lambda wc: work_racon_dir(wc.isolate)
    benchmark:
        bench("{isolate}", "racon_seed")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        cp "{input}" "{output}"
        """

rule racon_map:
    input:
        asm=lambda wc: stage_racon(wc.isolate, int(wc.i) - 1) if int(wc.i) > 1 else stage_racon0(wc.isolate),
        ont=filt_ont("{isolate}")
    output:
        paf=os.path.join(work_racon_dir("{isolate}"), "ont_vs_racon{i}.paf")
    params:
        polish=POLISH_MODE, rounds=RACON_ROUNDS,
        outdir=lambda wc: work_racon_dir(wc.isolate),
        mode=lambda wc: RUN_MODE[wc.isolate]
    log:
        os.path.join(logs_dir("{isolate}"), "racon_map_{i}.log")
    benchmark:
        bench("{isolate}", "racon_map_{i}")
    conda:
        "workflow/envs/assembly.yml"
    resources:
        mem_mb=POLISH_MEM_MB
    threads: THREADS
    message:
        "Racon map round {wildcards.i}"
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}"
        if [ "{params.mode}" = "illumina_only" ] || [ "{POLISH_MODE}" = "none" ]; then
          : > "{output.paf}"
        else
          minimap2 -t {threads} -x map-ont "{input.asm}" "{input.ont}" > "{output.paf}"
        fi
        ) > "{log}" 2>&1
        """

rule racon_consensus:
    input:
        ont=filt_ont("{isolate}"),
        paf=os.path.join(work_racon_dir("{isolate}"), "ont_vs_racon{i}.paf"),
        asm=lambda wc: stage_racon(wc.isolate, int(wc.i) - 1) if int(wc.i) > 1 else stage_racon0(wc.isolate)
    output:
        fasta=os.path.join(work_racon_dir("{isolate}"), "stage_racon{i}.fasta")
    params:
        polish=POLISH_MODE, rounds=RACON_ROUNDS,
        outdir=lambda wc: work_racon_dir(wc.isolate),
        mode=lambda wc: RUN_MODE[wc.isolate]
    log:
        os.path.join(logs_dir("{isolate}"), "racon_consensus_{i}.log")
    benchmark:
        bench("{isolate}", "racon_consensus_{i}")
    conda:
        "workflow/envs/assembly.yml"
    resources:
        mem_mb=POLISH_MEM_MB
    threads: THREADS
    message:
        "Racon consensus round {wildcards.i}"
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}"
        if [ "{params.mode}" = "illumina_only" ] || [ "{POLISH_MODE}" = "none" ]; then
          cp "{input.asm}" "{output.fasta}"
        else
          test -s "{input.paf}" || {{ echo "No alignments for requested Racon stage" >&2; exit 1; }}
          racon -t {threads} "{input.ont}" "{input.paf}" "{input.asm}" > "{output.fasta}"
        fi
        ) > "{log}" 2>&1
        """

rule medaka:
    input:
        prev=lambda wc: stage_racon(wc.isolate, RACON_ROUNDS) if RACON_ROUNDS > 0 else stage_unicycler(wc.isolate),
        ont=filt_ont("{isolate}")
    output:
        stage_medaka("{isolate}")
    params:
        polish=POLISH_MODE, rounds=RACON_ROUNDS,
        outdir=lambda wc: work_medaka_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate),
        mode=lambda wc: RUN_MODE[wc.isolate],
        model=lambda wc: ISO_OVERRIDES[wc.isolate]["medaka_model"]
    log:
        os.path.join(logs_dir("{isolate}"), "medaka.log")
    benchmark:
        bench("{isolate}", "medaka")
    conda:
        "workflow/envs/medaka.yml"
    resources:
        mem_mb=POLISH_MEM_MB
    threads: THREADS
    message:
        logmsg(13, "Medaka polishing")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"
        if [ "{params.mode}" = "illumina_only" ] || [ "{POLISH_MODE}" = "none" ] || [ "{POLISH_MODE}" = "racon" ]; then
          cp "{input.prev}" "{output}"
          exit 0
        fi
        if [ -n "{params.model}" ]; then
          medaka_consensus -t {threads} -i "{input.ont}" -d "{input.prev}" -o "{params.outdir}" -m "{params.model}"
        else
          medaka_consensus -t {threads} -i "{input.ont}" -d "{input.prev}" -o "{params.outdir}"
        fi
        ) > "{log}" 2>&1
        """

rule bwa_index:
    input:
        lambda wc: polishing_source(wc.isolate)
    output:
        ref=bwa_ref("{isolate}"), fai=bwa_fai("{isolate}"),
        amb=bwa_idx("{isolate}", "amb"), ann=bwa_idx("{isolate}", "ann"),
        bwt=bwa_idx("{isolate}", "bwt"), pac=bwa_idx("{isolate}", "pac"), sa=bwa_idx("{isolate}", "sa")
    params:
        outdir=lambda wc: work_mapping_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "bwa_index.log")
    benchmark: bench("{isolate}", "bwa_index")
    conda: "workflow/envs/polishing.yml"
    resources: mem_mb=POLISH_MEM_MB
    threads: 1
    message: "Index assembly for short-read polishing: {wildcards.isolate}"
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        (
          cp "{input}" "{output.ref}"
          samtools faidx "{output.ref}"
          bwa index "{output.ref}"
        ) > "{log}" 2>&1
        """

rule polypolish_align:
    input:
        ref=bwa_ref("{isolate}"), amb=bwa_idx("{isolate}", "amb"),
        ann=bwa_idx("{isolate}", "ann"), bwt=bwa_idx("{isolate}", "bwt"),
        pac=bwa_idx("{isolate}", "pac"), sa=bwa_idx("{isolate}", "sa"),
        r1=trim_r1("{isolate}"), r2=trim_r2("{isolate}")
    output:
        r1=polishing_sam("{isolate}", 1), r2=polishing_sam("{isolate}", 2)
    log: os.path.join(logs_dir("{isolate}"), "polypolish_align.log")
    benchmark: bench("{isolate}", "polypolish_align")
    conda: "workflow/envs/polishing.yml"
    resources: mem_mb=POLISH_MEM_MB, disk_mb=rule_memory("mapping_disk", 12000)
    threads: THREADS
    message: "Prepare separate all-alignment SAMs: {wildcards.isolate}"
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
          bwa mem -t {threads} -a "{input.ref}" "{input.r1}" > "{output.r1}"
          bwa mem -t {threads} -a "{input.ref}" "{input.r2}" > "{output.r2}"
        ) > "{log}" 2>&1
        """

rule illumina_map:
    input:
        ref=bwa_ref("{isolate}"), fai=bwa_fai("{isolate}"), amb=bwa_idx("{isolate}", "amb"),
        ann=bwa_idx("{isolate}", "ann"), bwt=bwa_idx("{isolate}", "bwt"),
        pac=bwa_idx("{isolate}", "pac"), sa=bwa_idx("{isolate}", "sa"),
        r1=trim_r1("{isolate}"), r2=trim_r2("{isolate}")
    output:
        bam=coord_bam("{isolate}"), bai=coord_bai("{isolate}"), qn=namesort_bam("{isolate}"),
        count_file=alignment_count("{isolate}")
    params:
        outdir=lambda wc: work_mapping_dir(wc.isolate),
        sort_extra=lambda wc, threads: max(0, threads - 1),
        sort_memory=lambda wc, threads: max(64, int(POLISH_MEM_MB * 0.7 / max(1, threads)))
    log: os.path.join(logs_dir("{isolate}"), "illumina_map.log")
    benchmark: bench("{isolate}", "illumina_map")
    conda: "workflow/envs/polishing.yml"
    resources: mem_mb=POLISH_MEM_MB, disk_mb=rule_memory("mapping_disk", 12000)
    threads: THREADS
    message: "Create paired-end mapping QC: {wildcards.isolate}"
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log}")" "$(dirname "{output.count_file}")"
        (
          tmp="{params.outdir}/.paired_qc.sam"
          trap 'rm -f "$tmp"' EXIT
          bwa mem -t {threads} -a -x intractg "{input.ref}" "{input.r1}" "{input.r2}" > "$tmp"
          samtools sort -@ {params.sort_extra} -m {params.sort_memory}M -O BAM -o "{output.bam}" "$tmp"
          samtools index "{output.bam}"
          samtools sort -n -@ {params.sort_extra} -m {params.sort_memory}M -O BAM "{output.bam}" -o "{output.qn}"
          samtools view -c "{output.qn}" > "{output.count_file}"
        ) > "{log}" 2>&1
        """

rule polypolish:
    input:
        ref=bwa_ref("{isolate}"), r1=polishing_sam("{isolate}", 1), r2=polishing_sam("{isolate}", 2)
    output:
        fasta=stage_polypolish("{isolate}"), status=polish_status("{isolate}")
    params:
        outdir=lambda wc: work_polypolish_dir(wc.isolate)
    log: os.path.join(logs_dir("{isolate}"), "polypolish.log")
    benchmark: bench("{isolate}", "polypolish")
    conda: "workflow/envs/polishing.yml"
    resources: mem_mb=POLISH_MEM_MB
    threads: 1
    message: "Polypolish: {wildcards.isolate}"
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log}")" "$(dirname "{output.status}")"
        (
          count1=$(samtools view -c -F 4 "{input.r1}")
          count2=$(samtools view -c -F 4 "{input.r2}")
          if [ "$((count1 + count2))" -eq 0 ]; then
            cp "{input.ref}" "{output.fasta}"
            printf '{{"status":"skipped","reason":"no_mapped_alignments","mapped_alignments":0}}\n' > "{output.status}"
          else
            tmp="{output.fasta}.tmp"
            trap 'rm -f "$tmp"' EXIT
            polypolish polish "{input.ref}" "{input.r1}" "{input.r2}" > "$tmp"
            test -s "$tmp"
            mv "$tmp" "{output.fasta}"
            printf '{{"status":"completed","reason":"","mapped_alignments":%s}}\n' "$((count1 + count2))" > "{output.status}"
          fi
        ) > "{log}" 2>&1
        """

rule polish_execution:
    input:
        assemblies=lambda wc: [stage_path(wc.isolate, x) for x in [base_stage(wc.isolate)] + planned_stages(wc.isolate)],
        polypolish_status=lambda wc: [polish_status(wc.isolate)] if needs_polypolish(wc.isolate) else []
    output: polish_execution("{isolate}")
    params:
        requested=POLISH_MODE, rounds=RACON_ROUNDS, polypolish=RUN_POLYPOLISH,
        mode=lambda wc: RUN_MODE[wc.isolate],
        source=lambda wc: final_source_for(wc.isolate),
        helper_sha=sha256_file(SUPPORT_FILE)
    message: "Record executed polishing stages: {wildcards.isolate}"
    benchmark:
        bench("{isolate}", "polish_execution")
    run:
        iso = wildcards.isolate
        planned = planned_stages(iso)
        executed = list(planned)
        skipped = {}
        if "polypolish" not in planned:
            skipped["polypolish"] = "disabled_or_incompatible_with_requested_mode"
        for stage in ("racon", "medaka"):
            if not any(x.startswith(stage) for x in planned):
                skipped[stage] = "disabled_or_incompatible_with_requested_mode"
        pp = {}
        if input.polypolish_status:
            pp = json.loads(Path(input.polypolish_status[0]).read_text())
            if pp["status"] == "skipped":
                executed.remove("polypolish")
                skipped["polypolish"] = pp["reason"]
        stages = [{"stage": name, "path": path, **require_fasta(path)}
                  for name, path in zip([base_stage(iso)] + planned, input.assemblies)]
        final_stage = executed[-1] if executed else base_stage(iso)
        write_json(output[0], {"schema": "hybas-polishing-v1", "workflow_version": VERSION,
            "requested_polish_mode": POLISH_MODE, "run_mode": RUN_MODE[iso],
            "planned_stages": planned, "executed_stages": executed, "skipped_stages": skipped,
            "final_source_stage": final_stage, "final_source_path": stage_path(iso, final_stage),
            "promoted_artifact_path": params.source,
            "final_source_sha256": sha256_file(params.source), "stages": stages,
            "polypolish": pp})

###############################################################################
# SECTION 18 — FINAL ASSEMBLY
###############################################################################

def final_source_for(i):
    stages = planned_stages(i)
    return stage_path(i, stages[-1] if stages else base_stage(i))

rule final_copy:
    input:
        source=lambda wc: final_source_for(wc.isolate), execution=polish_execution("{isolate}")
    output:
        report(final_asm("{isolate}"), caption="workflow/report/final_assembly.rst", category="Assembly", subcategory="Final")
    params:
        outdir=lambda wc: work_final_dir(wc.isolate)
    log: os.path.join(logs_dir("{isolate}"), "final_copy.log")
    message: "Promote validated final assembly: {wildcards.isolate}"
    benchmark:
        bench("{isolate}", "final_copy")
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        cp "{input.source}" "{output}.tmp"
        mv "{output}.tmp" "{output}"
        printf '%s\n' "Final assembly copied from {input.source}" > "{log}"
        """

###############################################################################
# SECTION 19 — QC / ANALYTICS OUTPUTS
###############################################################################

rule quast:
    input:
        final_asm("{isolate}")
    output:
        txt=report(quast_txt("{isolate}"), caption="workflow/report/quast_txt.rst", category="QC", subcategory="QUAST"),
        tsv=report(quast_tsv("{isolate}"), caption="workflow/report/quast_tsv.rst", category="QC", subcategory="QUAST")
    params:
        outdir=lambda wc: quast_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "quast.log")
    benchmark:
        bench("{isolate}", "quast")
    conda:
        "workflow/envs/qc.yml"
    threads: THREADS
    message:
        logmsg(18, "QUAST")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}"
        quast.py -o "{params.outdir}" -t {threads} "{input}"
        ) > "{log}" 2>&1
        """

rule busco:
    input:
        final_asm("{isolate}")
    output:
        txt=report(busco_txt("{isolate}"), caption="workflow/report/busco_txt.rst", category="QC", subcategory="BUSCO"),
        html=report(busco_html("{isolate}"), caption="workflow/report/busco_html.rst", category="QC", subcategory="BUSCO")
    params:
        enabled=RUN_BUSCO,
        outdir=lambda wc: busco_dir(wc.isolate),
        rundir=lambda wc: busco_run_dir(wc.isolate),
        logdir=lambda wc: logs_dir(wc.isolate),
        lineage=lambda wc: ISO_OVERRIDES[wc.isolate]["busco_lineage"]
    log:
        os.path.join(logs_dir("{isolate}"), "busco.log")
    benchmark:
        bench("{isolate}", "busco")
    conda:
        "workflow/envs/qc.yml"
    resources:
        mem_mb=rule_memory("busco", 4096)
    threads: THREADS
    message:
        logmsg(19, "BUSCO")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.logdir}"
        if [ "{RUN_BUSCO}" != "True" ]; then
          mkdir -p "{params.rundir}"
          : > "{output.txt}"
          : > "{output.html}"
          exit 0
        fi
        busco -i "{input}" \
          -o "busco_{wildcards.isolate}" \
          --out_path "{params.outdir}" \
          -m genome \
          -l "{BUSCO_DB_DIR}/{params.lineage}" --offline \
          --download_path "{BUSCO_DB_DIR}" \
          -c {threads} \
          --force
        sfile=$(ls -1 "{params.rundir}"/short_summary*.txt 2>/dev/null | head -n1 || true)
        hfile=$(ls -1 "{params.rundir}"/short_summary*.html 2>/dev/null | head -n1 || true)
        if [ -n "$sfile" ] && [ -s "$sfile" ]; then if [ "$sfile" != "{output.txt}" ]; then cp -f "$sfile" "{output.txt}"; fi; else echo "BUSCO completed without a nonempty short summary" >&2; exit 1; fi
        if [ -n "$hfile" ] && [ -s "$hfile" ]; then if [ "$hfile" != "{output.html}" ]; then cp -f "$hfile" "{output.html}"; fi; else printf '<html><body>BUSCO does not provide an HTML summary in this version. See short_summary.txt.</body></html>\n' > "{output.html}"; fi
        ) > "{log}" 2>&1
        """

rule seqkit_stats:
    input:
        raw_r1=merged_r1("{isolate}"),
        raw_r2=merged_r2("{isolate}"),
        clean_r1=clean_r1("{isolate}"),
        clean_r2=clean_r2("{isolate}"),
        trim_r1=trim_r1("{isolate}"),
        trim_r2=trim_r2("{isolate}"),
        ont_raw=merged_ont("{isolate}"),
        ont_clean=clean_ont_path("{isolate}"),
        ont_filt=filt_ont("{isolate}"),
        asm=final_asm("{isolate}")
    output:
        illumina=report(
            seqkit_illumina("{isolate}"),
            caption="workflow/report/seqkit_tsv.rst",
            category="QC",
            subcategory="seqkit",
        ),
        ont=report(
            seqkit_ont("{isolate}"),
            caption="workflow/report/seqkit_tsv.rst",
            category="QC",
            subcategory="seqkit",
        ),
        asm=report(
            seqkit_assembly("{isolate}"),
            caption="workflow/report/seqkit_tsv.rst",
            category="QC",
            subcategory="seqkit",
        )
    params:
        outdir=lambda wc: seqkit_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "seqkit_stats.log")
    conda:
        "workflow/envs/qc.yml"
    threads: 1
    message:
        logmsg(20, "seqkit summaries")
    benchmark:
        bench("{isolate}", "seqkit_stats")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}"

        {{
          echo -e "stage\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len"
          if [ -s "{input.raw_r1}" ];   then seqkit stats -a -T "{input.raw_r1}"   | awk 'NR==2{{print "raw_R1\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}'; fi
          if [ -s "{input.raw_r2}" ];   then seqkit stats -a -T "{input.raw_r2}"   | awk 'NR==2{{print "raw_R2\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}'; fi
          if [ -s "{input.clean_r1}" ]; then seqkit stats -a -T "{input.clean_r1}" | awk 'NR==2{{print "clean_R1\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}'; fi
          if [ -s "{input.clean_r2}" ]; then seqkit stats -a -T "{input.clean_r2}" | awk 'NR==2{{print "clean_R2\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}'; fi
          if [ -s "{input.trim_r1}" ];  then seqkit stats -a -T "{input.trim_r1}"  | awk 'NR==2{{print "trim_R1\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}'; fi
          if [ -s "{input.trim_r2}" ];  then seqkit stats -a -T "{input.trim_r2}"  | awk 'NR==2{{print "trim_R2\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}'; fi
        }} > "{output.illumina}"

        {{
          echo -e "stage\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len"
          if [ -s "{input.ont_raw}" ];   then seqkit stats -a -T "{input.ont_raw}"   | awk 'NR==2{{print "raw_ONT\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}'; fi
          if [ -s "{input.ont_clean}" ]; then seqkit stats -a -T "{input.ont_clean}" | awk 'NR==2{{print "clean_ONT\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}'; fi
          if [ -s "{input.ont_filt}" ];  then seqkit stats -a -T "{input.ont_filt}"  | awk 'NR==2{{print "filtlong_ONT\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}'; fi
        }} > "{output.ont}"

        {{
          echo -e "stage\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len"
          seqkit stats -a -T "{input.asm}" | awk 'NR==2{{print "final_assembly\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}'
        }} > "{output.asm}"

        ) > "{log}" 2>&1
        """

rule coverage_estimate:
    input:
        ont=filt_ont("{isolate}"),
        asm=final_asm("{isolate}")
    output:
        note=report(
            coverage_note("{isolate}"),
            caption="workflow/report/coverage_txt.rst",
            category="QC",
            subcategory="Coverage",
        ),
        flag=coverage_checked("{isolate}")
    params:
        outdir=lambda wc: coverage_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "coverage_estimate.log")
    conda:
        "workflow/envs/qc.yml"
    threads: 1
    message:
        logmsg(21, "ONT coverage estimate")
    benchmark:
        bench("{isolate}", "coverage_estimate")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}"
        if [ -s "{input.ont}" ]; then
          ont_bp=$(seqkit stats -a -T "{input.ont}" | awk 'NR==2 {{gsub(/,/,"",$5); print $5}}')
        else
          ont_bp=0
        fi
        asm_bp=$(seqkit stats -a -T "{input.asm}" | awk 'NR==2 {{gsub(/,/,"",$5); print $5}}')
        cov=$(awk -v a="$ont_bp" -v b="$asm_bp" 'BEGIN{{if(b>0) printf "%.2f", a/b; else print 0}}')
        printf "field\tvalue\nont_bp\t%s\nassembly_bp\t%s\nont_coverage_x\t%s\n" "$ont_bp" "$asm_bp" "$cov" > "{output.note}"
        : > "{output.flag}"
        ) > "{log}" 2>&1
        """

rule multiqc:
    input:
        custom=multiqc_custom_tsv("{isolate}"),
        config=os.path.join(SOURCE_ROOT, "workflow", "config", "multiqc_config.yaml"),
        quast_txt=quast_txt("{isolate}"),
        quast_tsv=quast_tsv("{isolate}"),
        busco=busco_txt("{isolate}"),
        fastp_html=os.path.join(fastp_dir("{isolate}"), "fastp.html"),
        fastp_json=os.path.join(fastp_dir("{isolate}"), "fastp.json"),
        k1=os.path.join(kraken_dir("{isolate}"), "illumina.report.txt"),
        k2=os.path.join(kraken_dir("{isolate}"), "nanopore.report.txt"),
        s1=seqkit_illumina("{isolate}"),
        s2=seqkit_ont("{isolate}"),
        s3=seqkit_assembly("{isolate}"),
        cov=coverage_note("{isolate}"),
        qc_assess=os.path.join(metadata_dir("{isolate}"), "qc_assessment.tsv"),
        fastqc_dir=fastqc_pre_dir("{isolate}"),
        nanoplot_raw=nanoplot_raw_dir("{isolate}"),
        nanoplot_filt=nanoplot_filt_dir("{isolate}")
    output:
        report(
            multiqc_html("{isolate}"),
            caption="workflow/report/multiqc.rst",
            category="QC",
            subcategory="MultiQC",
        )
    params:
        enabled=RUN_MULTIQC,
        outdir=lambda wc: multiqc_dir(wc.isolate),
        isodir=lambda wc: iso_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "multiqc.log")
    benchmark:
        bench("{isolate}", "multiqc")
    conda:
        "workflow/envs/qc.yml"
    threads: 1
    message:
        logmsg(22, "MultiQC")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}"
        if [ "{RUN_MULTIQC}" != "True" ]; then
          : > "{output}"
          exit 0
        fi
        multiqc -f -c "{input.config}" -n multiqc_report.html -o "{params.outdir}" "{params.isodir}/reports" "{params.isodir}/metadata" --ignore "*/multiqc/*" --ignore "*/safety/*"
        ) > "{log}" 2>&1
        """

###############################################################################
# SECTION 20 — ANNOTATION
###############################################################################

rule prokka:
    input:
        final_asm("{isolate}")
    output:
        gff=report(prokka_gff("{isolate}"), caption="workflow/report/prokka_gff.rst", category="Annotation", subcategory="Prokka"),
        faa=report(prokka_faa("{isolate}"), caption="workflow/report/prokka_faa.rst", category="Annotation", subcategory="Prokka")
    params:
        enabled=RUN_PROKKA,
        outdir=lambda wc: annotation_dir(wc.isolate)
    log:
        os.path.join(logs_dir("{isolate}"), "prokka.log")
    benchmark:
        bench("{isolate}", "prokka")
    conda:
        "workflow/envs/annotation.yml" if RUN_PROKKA else "workflow/envs/controller.yml"
    resources:
        mem_mb=rule_memory("prokka", 4096)
    threads: THREADS
    message:
        logmsg(23, "Prokka")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        mkdir -p "{params.outdir}"
        if [ "{RUN_PROKKA}" != "True" ]; then
          : > "{output.gff}"
          : > "{output.faa}"
          exit 0
        fi
        prokka "{input}" --outdir "{params.outdir}" --prefix "{wildcards.isolate}" --cpus {threads} --force
        ) > "{log}" 2>&1
        """

###############################################################################
# SECTION 21 — TOOL VERSION SNAPSHOT
###############################################################################

rule tool_versions:
    input:
        assembly_list=lambda wc: active_env_list("assembly"),
        library_list=lambda wc: active_env_list("library_assembly") if wc.isolate in LIBRARY_SAMPLES else [],
        qc_list=lambda wc: active_env_list("qc"),
        polishing_list=lambda wc: active_env_list("polishing"),
        medaka_list=lambda wc: active_env_list("medaka"),
        annotation_list=lambda wc: active_env_list("annotation")
    output:
        report(
            tool_versions_tsv("{isolate}"),
            caption="workflow/report/tool_versions_tsv.rst",
            category="Run metadata",
            subcategory="Tool versions",
        )
    benchmark:
        bench("{isolate}", "tool_versions")
    run:
        iso = wildcards.isolate
        os.makedirs(metadata_dir(iso), exist_ok=True)

        def parse_conda_list_version(path, names):
            if isinstance(names, str):
                names = [names]
            if isinstance(path, (list, tuple)):
                if not path: return "NA"
                path = path[0]
            if not os.path.exists(path):
                return "NA"

            wanted = {str(n).lower() for n in names}
            with open(path, "r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    s = line.strip()
                    if not s or s.startswith("#"):
                        continue
                    parts = s.split()
                    if len(parts) >= 2 and parts[0].lower() in wanted:
                        return parts[1]
            return "NA"

        def first_version(*candidates):
            for value in candidates:
                if value and value != "NA":
                    return value
            return "NA"

        versions = {
            "snakemake": __import__("importlib.metadata", fromlist=["version"]).version("snakemake"),

            "unicycler": parse_conda_list_version(input.assembly_list, "unicycler"),
            "spades": parse_conda_list_version(input.library_list if iso in LIBRARY_SAMPLES else input.assembly_list, "spades"),
            "racon": parse_conda_list_version(input.assembly_list, "racon"),
            "minimap2": parse_conda_list_version(input.assembly_list, "minimap2"),
            "filtlong": parse_conda_list_version(input.assembly_list, "filtlong"),
            "pigz": parse_conda_list_version(input.assembly_list, "pigz"),
            "fastp": first_version(
                parse_conda_list_version(input.assembly_list, "fastp"),
                parse_conda_list_version(input.qc_list, "fastp"),
            ),

            "fastqc": parse_conda_list_version(input.qc_list, "fastqc"),
            "kraken2": parse_conda_list_version(input.qc_list, "kraken2"),
            "krakentools": parse_conda_list_version(input.qc_list, ["krakentools", "krakentools.py", "krakentools-py"]),
            "busco": parse_conda_list_version(input.qc_list, "busco"),
            "quast": parse_conda_list_version(input.qc_list, "quast"),
            "multiqc": parse_conda_list_version(input.qc_list, "multiqc"),
            "nanoplot": parse_conda_list_version(input.qc_list, ["nanoplot", "NanoPlot"]),
            "seqkit": parse_conda_list_version(input.qc_list, "seqkit"),

            "bwa": parse_conda_list_version(input.polishing_list, "bwa"),
            "samtools": parse_conda_list_version(input.polishing_list, "samtools"),
            "polypolish": parse_conda_list_version(input.polishing_list, "polypolish"),

            "medaka": parse_conda_list_version(input.medaka_list, "medaka"),

            "prokka": parse_conda_list_version(input.annotation_list, "prokka"),
            "blast": parse_conda_list_version(input.annotation_list, ["blast", "blast-legacy"]),
            "barrnap": parse_conda_list_version(input.annotation_list, "barrnap"),
        }

        with open(output[0], "w", encoding="utf-8") as out:
            out.write("tool\tversion\n")
            for tool, version in versions.items():
                out.write(f"{tool}\t{version}\n")

###############################################################################
# SECTION 22 — CHECKSUMS
###############################################################################

rule checksums:
    input:
        final=final_asm("{isolate}"),
        runmode=run_mode_file("{isolate}"),
        summary=run_summary_json("{isolate}")
    output:
        report(
            checksums_tsv("{isolate}"),
            caption="workflow/report/checksums_tsv.rst",
            category="Run metadata",
            subcategory="Checksums",
        )
    benchmark:
        bench("{isolate}", "checksums")
    run:
        os.makedirs(metadata_dir(wildcards.isolate), exist_ok=True)

        files = []
        for p in illumina_r1_files(wildcards.isolate):
            files.append(("raw_illumina_r1", p))
        for p in illumina_r2_files(wildcards.isolate):
            files.append(("raw_illumina_r2", p))
        for p in nanopore_files(wildcards.isolate):
            files.append(("raw_ont", p))

        files.extend([
            ("merged_r1", merged_r1(wildcards.isolate)),
            ("merged_r2", merged_r2(wildcards.isolate)),
            ("merged_ont", merged_ont(wildcards.isolate)),
            ("final_assembly", final_asm(wildcards.isolate)),
            ("run_mode", run_mode_file(wildcards.isolate)),
            ("run_summary_json", run_summary_json(wildcards.isolate)),
        ])

        with open(output[0], "w", encoding="utf-8") as out:
            out.write("label\tpath\tsha256\n")
            for label, path in files:
                if os.path.exists(path):
                    out.write(f"{label}\t{path}\t{sha256_file(path)}\n")
                else:
                    out.write(f"{label}\t{path}\tNA\n")

###############################################################################
# SECTION 23 — RUN SUMMARY
###############################################################################

rule run_summary:
    params:
        version=VERSION, config=EFFECTIVE_CONFIG, helper_sha=sha256_file(SUPPORT_FILE)
    input:
        final=final_asm("{isolate}"),
        mode=run_mode_file("{isolate}"),
        seq_illumina=seqkit_illumina("{isolate}"),
        seq_ont=seqkit_ont("{isolate}"),
        seq_asm=seqkit_assembly("{isolate}"),
        quast=quast_tsv("{isolate}"),
        busco=busco_txt("{isolate}"),
        cov=coverage_note("{isolate}"),
        skip=skip_reasons_tsv("{isolate}"),
        bam=lambda wc: [namesort_bam(wc.isolate)] if needs_polypolish(wc.isolate) else [],
        mapping_count=lambda wc: [alignment_count(wc.isolate)] if needs_polypolish(wc.isolate) else [],
        execution=polish_execution("{isolate}"),
        kr_illumina=os.path.join(kraken_dir("{isolate}"), "illumina.report.txt"),
        kr_ont=os.path.join(kraken_dir("{isolate}"), "nanopore.report.txt"),
        kraken_decision=kraken_decision_tsv("{isolate}"),
        fastp_html=os.path.join(fastp_dir("{isolate}"), "fastp.html"),
        fastp_json=os.path.join(fastp_dir("{isolate}"), "fastp.json"),
        tools=tool_versions_tsv("{isolate}")
    output:
        tsv=report(run_summary_tsv("{isolate}"), caption="workflow/report/run_summary_tsv.rst", category="Run metadata", subcategory="TSV"),
        json=report(run_summary_json("{isolate}"), caption="workflow/report/run_summary_json.rst", category="Run metadata", subcategory="JSON")
    benchmark:
        bench("{isolate}", "run_summary")
    run:
        import subprocess

        iso = wildcards.isolate
        os.makedirs(metadata_dir(iso), exist_ok=True)

        def read_text(path):
            return open(path, "r", encoding="utf-8", errors="ignore").read() if os.path.exists(path) else ""

        def parse_tsv_dict(path, key_col=0, val_col=1, has_header=True):
            data = {}
            if not os.path.exists(path):
                return data
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                lines = fh.read().splitlines()
            if has_header and lines:
                lines = lines[1:]
            for line in lines:
                parts = line.split("\t")
                if len(parts) > max(key_col, val_col):
                    data[parts[key_col]] = parts[val_col]
            return data

        def parse_seqkit_table(path):
            rows = []
            if not os.path.exists(path):
                return rows
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                lines = [x.rstrip("\n") for x in fh if x.strip()]
            if len(lines) < 2:
                return rows
            header = lines[0].split("\t")
            for line in lines[1:]:
                vals = line.split("\t")
                row = {}
                for i, h in enumerate(header):
                    row[h] = vals[i] if i < len(vals) else ""
                rows.append(row)
            return rows

        def seqkit_stage_value(path, stage, field):
            for row in parse_seqkit_table(path):
                if row.get("stage") == stage:
                    return row.get(field, "NA").replace(",", "")
            return "NA"

        def safe_fraction_local(num, den, ndigits=4):
            try:
                num = float(str(num).replace(",", ""))
                den = float(str(den).replace(",", ""))
                if den > 0:
                    return f"{num/den:.{ndigits}f}"
            except Exception:
                pass
            return "NA"

        def contig_count(path):
            if not os.path.exists(path):
                return "NA"
            c = 0
            with open(path, "r", encoding="utf-8", errors="ignore") as h:
                for line in h:
                    if line.startswith(">"):
                        c += 1
            return str(c)

        def quast_metric(path, metric_name):
            return read_quast_metric(path, metric_name)

        def parse_busco(path):
            txt = read_text(path)
            out = {
                "summary_line": "NA",
                "complete_percent": "NA",
                "single_copy_percent": "NA",
                "duplicated_percent": "NA",
                "fragmented_percent": "NA",
                "missing_percent": "NA",
                "n_buscos": "NA",
            }

            m = re.search(
                r"C:(?P<C>[0-9.]+)%\[S:(?P<S>[0-9.]+)%,D:(?P<D>[0-9.]+)%\],F:(?P<F>[0-9.]+)%,M:(?P<M>[0-9.]+)%,n:(?P<n>[0-9]+)",
                txt
            )
            if m:
                out["summary_line"] = m.group(0)
                out["complete_percent"] = m.group("C")
                out["single_copy_percent"] = m.group("S")
                out["duplicated_percent"] = m.group("D")
                out["fragmented_percent"] = m.group("F")
                out["missing_percent"] = m.group("M")
                out["n_buscos"] = m.group("n")
            else:
                for line in txt.splitlines():
                    if "C:" in line:
                        out["summary_line"] = line.strip()
                        break
            return out

        def parse_cov(path):
            return parse_tsv_dict(path)

        def parse_tool_version(path, tool_name):
            if not os.path.exists(path):
                return "NA"
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                lines = fh.read().splitlines()
            for line in lines[1:]:
                parts = line.split("\t")
                if len(parts) >= 2 and parts[0] == tool_name:
                    return parts[1]
            return "NA"

        illumina_raw_r1_reads = seqkit_stage_value(input.seq_illumina, "raw_R1", "num_seqs")
        illumina_raw_r2_reads = seqkit_stage_value(input.seq_illumina, "raw_R2", "num_seqs")
        illumina_clean_r1_reads = seqkit_stage_value(input.seq_illumina, "clean_R1", "num_seqs")
        illumina_clean_r2_reads = seqkit_stage_value(input.seq_illumina, "clean_R2", "num_seqs")
        illumina_trim_r1_reads = seqkit_stage_value(input.seq_illumina, "trim_R1", "num_seqs")
        illumina_trim_r2_reads = seqkit_stage_value(input.seq_illumina, "trim_R2", "num_seqs")

        illumina_raw_r1_bp = seqkit_stage_value(input.seq_illumina, "raw_R1", "sum_len")
        illumina_raw_r2_bp = seqkit_stage_value(input.seq_illumina, "raw_R2", "sum_len")
        illumina_clean_r1_bp = seqkit_stage_value(input.seq_illumina, "clean_R1", "sum_len")
        illumina_clean_r2_bp = seqkit_stage_value(input.seq_illumina, "clean_R2", "sum_len")
        illumina_trim_r1_bp = seqkit_stage_value(input.seq_illumina, "trim_R1", "sum_len")
        illumina_trim_r2_bp = seqkit_stage_value(input.seq_illumina, "trim_R2", "sum_len")

        ont_raw_reads = seqkit_stage_value(input.seq_ont, "raw_ONT", "num_seqs")
        ont_clean_reads = seqkit_stage_value(input.seq_ont, "clean_ONT", "num_seqs")
        ont_filt_reads = seqkit_stage_value(input.seq_ont, "filtlong_ONT", "num_seqs")

        ont_raw_bp = seqkit_stage_value(input.seq_ont, "raw_ONT", "sum_len")
        ont_clean_bp = seqkit_stage_value(input.seq_ont, "clean_ONT", "sum_len")
        ont_filt_bp = seqkit_stage_value(input.seq_ont, "filtlong_ONT", "sum_len")

        assembly_bp = seqkit_stage_value(input.seq_asm, "final_assembly", "sum_len")
        assembly_num_seqs = seqkit_stage_value(input.seq_asm, "final_assembly", "num_seqs")
        assembly_min_len = seqkit_stage_value(input.seq_asm, "final_assembly", "min_len")
        assembly_avg_len = seqkit_stage_value(input.seq_asm, "final_assembly", "avg_len")
        assembly_max_len = seqkit_stage_value(input.seq_asm, "final_assembly", "max_len")

        busco_stats = parse_busco(input.busco)
        cov_stats = parse_cov(input.cov)
        kraken_decision_stats = parse_tsv_dict(input.kraken_decision)

        execution = json.loads(Path(input.execution).read_text())
        polypolish_alignments = execution.get("polypolish", {}).get("mapped_alignments", "NA")
        mapping_alignment_count = Path(input.mapping_count[0]).read_text().strip() if input.mapping_count else "NA"

        summary = {
            "isolate": iso,
            "run_mode": RUN_MODE[iso],
            "workflow_version": WORKFLOW_VERSION,
            "polish_mode_requested": POLISH_MODE,
            "requested_polish_mode": POLISH_MODE,
            "executed_stages": execution["executed_stages"],
            "skipped_stages": execution["skipped_stages"],
            "final_source_stage": execution["final_source_stage"],
            "final_source_sha256": execution["final_source_sha256"],
            "mapping_alignment_count": mapping_alignment_count,

            "run_kraken": KRAKEN_PROFILE_ENABLED,
            "kraken_decision_mode": kraken_decision_stats.get("kraken_decision_mode", "NA"),
            "kraken_decision": kraken_decision_stats.get("kraken_decision", "NA"),
            "kraken_decision_reason": kraken_decision_stats.get("kraken_decision_reason", "NA"),

            "kraken_manual_target_taxid": kraken_decision_stats.get("manual_target_taxid", "NA"),

            "kraken_selected_target_taxid": kraken_decision_stats.get("selected_target_taxid", "NA"),
            "kraken_selected_target_name": kraken_decision_stats.get("selected_target_name", "NA"),
            "kraken_selected_target_rank": kraken_decision_stats.get("selected_target_rank", "NA"),
            "kraken_selected_target_source": kraken_decision_stats.get("selected_target_source", "NA"),
            "kraken_selected_target_percent_for_decision": kraken_decision_stats.get("selected_target_percent_for_decision", "NA"),
            "kraken_illumina_selected_target_percent": kraken_decision_stats.get("illumina_selected_target_percent", "NA"),
            "kraken_ont_selected_target_percent": kraken_decision_stats.get("ont_selected_target_percent", "NA"),

            "kraken_inferred_target_taxid": kraken_decision_stats.get("inferred_target_taxid", "NA"),
            "kraken_inferred_target_name": kraken_decision_stats.get("inferred_target_name", "NA"),
            "kraken_inferred_target_rank": kraken_decision_stats.get("inferred_target_rank", "NA"),
            "kraken_inferred_target_percent": kraken_decision_stats.get("inferred_target_percent", "NA"),
            "kraken_target_selection_reason": kraken_decision_stats.get("target_selection_reason", "NA"),

            "kraken_runner_up_taxid": kraken_decision_stats.get("runner_up_taxid", "NA"),
            "kraken_runner_up_name": kraken_decision_stats.get("runner_up_name", "NA"),
            "kraken_runner_up_percent": kraken_decision_stats.get("runner_up_percent", "NA"),
            "kraken_selected_vs_runner_up_delta_percent": kraken_decision_stats.get("selected_vs_runner_up_delta_percent", "NA"),
            "kraken_platform_agreement_score": kraken_decision_stats.get("platform_agreement_score", "NA"),
            "kraken_platform_agreement_reason": kraken_decision_stats.get("platform_agreement_reason", "NA"),
            "kraken_decision_confidence": kraken_decision_stats.get("decision_confidence", "NA"),
            "kraken_awareness_state": kraken_decision_stats.get("awareness_state", "NA"),

            "run_busco": RUN_BUSCO,
            "run_prokka": RUN_PROKKA,
            "run_multiqc": RUN_MULTIQC,
            "run_polypolish": RUN_POLYPOLISH,

            "target_taxid": ISO_OVERRIDES[iso]["target_taxid"],
            "busco_lineage": ISO_OVERRIDES[iso]["busco_lineage"],
            "medaka_model": ISO_OVERRIDES[iso]["medaka_model"],

            "input_illumina_r1_files": len(illumina_r1_files(iso)),
            "input_illumina_r2_files": len(illumina_r2_files(iso)),
            "input_ont_files": len(nanopore_files(iso)),

            "illumina_raw_r1_reads": illumina_raw_r1_reads,
            "illumina_raw_r2_reads": illumina_raw_r2_reads,
            "illumina_clean_r1_reads": illumina_clean_r1_reads,
            "illumina_clean_r2_reads": illumina_clean_r2_reads,
            "illumina_trim_r1_reads": illumina_trim_r1_reads,
            "illumina_trim_r2_reads": illumina_trim_r2_reads,

            "illumina_raw_r1_bp": illumina_raw_r1_bp,
            "illumina_raw_r2_bp": illumina_raw_r2_bp,
            "illumina_clean_r1_bp": illumina_clean_r1_bp,
            "illumina_clean_r2_bp": illumina_clean_r2_bp,
            "illumina_trim_r1_bp": illumina_trim_r1_bp,
            "illumina_trim_r2_bp": illumina_trim_r2_bp,

            "illumina_r1_kraken_retention_fraction": safe_fraction_local(illumina_clean_r1_bp, illumina_raw_r1_bp),
            "illumina_r2_kraken_retention_fraction": safe_fraction_local(illumina_clean_r2_bp, illumina_raw_r2_bp),
            "illumina_r1_trim_retention_fraction": safe_fraction_local(illumina_trim_r1_bp, illumina_clean_r1_bp),
            "illumina_r2_trim_retention_fraction": safe_fraction_local(illumina_trim_r2_bp, illumina_clean_r2_bp),

            "ont_raw_reads": ont_raw_reads,
            "ont_clean_reads": ont_clean_reads,
            "ont_filt_reads": ont_filt_reads,
            "ont_raw_bp": ont_raw_bp,
            "ont_clean_bp": ont_clean_bp,
            "ont_filt_bp": ont_filt_bp,
            "ont_kraken_retention_fraction": safe_fraction_local(ont_clean_bp, ont_raw_bp),
            "ont_filtlong_retention_fraction": safe_fraction_local(ont_filt_bp, ont_clean_bp),

            "assembly_bp": assembly_bp,
            "assembly_num_seqs_seqkit": assembly_num_seqs,
            "assembly_contig_count_fasta": contig_count(input.final),
            "assembly_min_len": assembly_min_len,
            "assembly_avg_len": assembly_avg_len,
            "assembly_max_len": assembly_max_len,

            "quast_total_length": quast_metric(input.quast, "Total length"),
            "quast_gc_percent": quast_metric(input.quast, "GC (%)"),
            "quast_n50": quast_metric(input.quast, "N50"),
            "quast_l50": quast_metric(input.quast, "L50"),
            "quast_num_contigs": quast_metric(input.quast, "# contigs"),
            "quast_largest_contig": quast_metric(input.quast, "Largest contig"),

            "busco_summary": busco_stats["summary_line"],
            "busco_complete_percent": busco_stats["complete_percent"],
            "busco_single_copy_percent": busco_stats["single_copy_percent"],
            "busco_duplicated_percent": busco_stats["duplicated_percent"],
            "busco_fragmented_percent": busco_stats["fragmented_percent"],
            "busco_missing_percent": busco_stats["missing_percent"],
            "busco_n_buscos": busco_stats["n_buscos"],

            "estimated_ont_coverage_x": cov_stats.get("ont_coverage_x", "NA"),
            "coverage_ont_bp": cov_stats.get("ont_bp", "NA"),
            "coverage_assembly_bp": cov_stats.get("assembly_bp", "NA"),

            "polypolish_alignment_count": polypolish_alignments,

            "kraken_illumina_report_present": os.path.exists(input.kr_illumina),
            "kraken_ont_report_present": os.path.exists(input.kr_ont),
            "assembler": base_stage(iso),
            "input_strategy": LIBRARY_SAMPLES[iso]["selection"]["selected_mode"] if iso in LIBRARY_SAMPLES else RUN_MODE[iso],
            "library_count": len(LIBRARY_SAMPLES[iso]["libraries"]) if iso in LIBRARY_SAMPLES else None,
            "tool_spades_version": parse_tool_version(input.tools, "spades"),
            "fastp_html_present": os.path.exists(input.fastp_html),
            "fastp_json_present": os.path.exists(input.fastp_json),

            "tool_snakemake_version": parse_tool_version(input.tools, "snakemake"),
            "tool_unicycler_version": parse_tool_version(input.tools, "unicycler"),
            "tool_racon_version": parse_tool_version(input.tools, "racon"),
            "tool_minimap2_version": parse_tool_version(input.tools, "minimap2"),
            "tool_filtlong_version": parse_tool_version(input.tools, "filtlong"),
            "tool_fastqc_version": parse_tool_version(input.tools, "fastqc"),
            "tool_kraken2_version": parse_tool_version(input.tools, "kraken2"),
            "tool_krakentools_version": parse_tool_version(input.tools, "krakentools"),
            "tool_fastp_version": parse_tool_version(input.tools, "fastp"),
            "tool_nanoplot_version": parse_tool_version(input.tools, "nanoplot"),
            "tool_seqkit_version": parse_tool_version(input.tools, "seqkit"),
            "tool_busco_version": parse_tool_version(input.tools, "busco"),
            "tool_quast_version": parse_tool_version(input.tools, "quast"),
            "tool_multiqc_version": parse_tool_version(input.tools, "multiqc"),
            "tool_bwa_version": parse_tool_version(input.tools, "bwa"),
            "tool_samtools_version": parse_tool_version(input.tools, "samtools"),
            "tool_polypolish_version": parse_tool_version(input.tools, "polypolish"),
            "tool_medaka_version": parse_tool_version(input.tools, "medaka"),
            "tool_prokka_version": parse_tool_version(input.tools, "prokka"),
            "tool_barrnap_version": parse_tool_version(input.tools, "barrnap"),

            "final_source": final_source_for(iso),
            "kraken_db_dir": KRAKEN_DB_DIR,
            "busco_lineage_dir": busco_lineage_dir(ISO_OVERRIDES[iso]["busco_lineage"]),
            "timestamp": datetime.now().isoformat(timespec="seconds"),
        }

        with open(output.tsv, "w", encoding="utf-8") as out:
            out.write("metric\tvalue\n")
            for k, v in summary.items():
                out.write(f"{k}\t{v}\n")

        with open(output.json, "w", encoding="utf-8") as out:
            json.dump(summary, out, indent=2, sort_keys=True)
            
###############################################################################
# SECTION 24 — QC ASSESSMENT AND MULTIQC CUSTOM CONTENT
###############################################################################

rule qc_assessment:
    params:
        busco=RUN_BUSCO, kraken=KRAKEN_PROFILE_ENABLED, mode=lambda wc: RUN_MODE[wc.isolate], helper_sha=sha256_file(SUPPORT_FILE)
    input:
        run_summary=run_summary_tsv("{isolate}"),
        quast=quast_tsv("{isolate}"),
        busco=busco_txt("{isolate}"),
        cov=coverage_note("{isolate}"),
        kr_illumina=os.path.join(kraken_dir("{isolate}"), "illumina.report.txt"),
        kr_ont=os.path.join(kraken_dir("{isolate}"), "nanopore.report.txt"),
        fastp_html=os.path.join(fastp_dir("{isolate}"), "fastp.html"),
        seq_illumina=seqkit_illumina("{isolate}"),
        seq_ont=seqkit_ont("{isolate}"),
        seq_asm=seqkit_assembly("{isolate}")
    output:
        report(
            os.path.join(metadata_dir("{isolate}"), "qc_assessment.tsv"),
            caption="workflow/report/qc_assessment_tsv.rst",
            category="Run metadata",
            subcategory="QC assessment",
        )
    benchmark:
        bench("{isolate}", "qc_assessment")
    run:
        iso = wildcards.isolate
        os.makedirs(metadata_dir(iso), exist_ok=True)

        def read_text(path):
            return open(path, "r", encoding="utf-8", errors="ignore").read() if os.path.exists(path) else ""

        def read_metric_table(path):
            data = {}
            for line in read_text(path).splitlines()[1:]:
                parts = line.split("\t", 1)
                if len(parts) == 2:
                    data[parts[0]] = parts[1]
            return data

        def busco_complete_percent(path):
            for line in read_text(path).splitlines():
                if "C:" in line:
                    try:
                        first = line.split("%")[0]
                        return first.replace("C:", "").strip()
                    except Exception:
                        return "NA"
            return "NA"

        def quast_metric(path, metric_name):
            return read_quast_metric(path, metric_name)

        runsum = read_metric_table(input.run_summary)
        cov = read_metric_table(input.cov)

        assembly_bp = runsum.get("assembly_bp", "NA")
        contigs = runsum.get("assembly_contig_count_fasta", "NA")
        ont_cov = cov.get("ont_coverage_x", "NA")
        busco_c = busco_complete_percent(input.busco)
        n50 = quast_metric(input.quast, "N50")

        rows = []

        rows.append(("assembly_present", "PASS" if os.path.exists(final_asm(iso)) else "FAIL", "final assembly file existence"))
        rows.append(("quast_present", "PASS" if os.path.exists(input.quast) else "FAIL", "QUAST report presence"))
        rows.append(("busco_present", ("PASS" if os.path.getsize(input.busco) > 0 else "FAIL") if RUN_BUSCO else "SKIP", "BUSCO short summary presence"))
        rows.append(("kraken_illumina_report_present", ("PASS" if os.path.getsize(input.kr_illumina) > 0 else "WARN") if KRAKEN_PROFILE_ENABLED and RUN_MODE[iso] != "ont_only" else "SKIP", "Kraken Illumina report presence"))
        rows.append(("kraken_nanopore_report_present", ("PASS" if os.path.getsize(input.kr_ont) > 0 else "WARN") if KRAKEN_PROFILE_ENABLED and RUN_MODE[iso] != "illumina_only" else "SKIP", "Kraken ONT report presence"))
        rows.append(("fastp_report_present", ("PASS" if os.path.getsize(input.fastp_html) > 0 else "WARN") if RUN_MODE[iso] != "ont_only" else "SKIP", "fastp HTML report presence"))
        rows.append(("seqkit_illumina_present", "PASS" if os.path.exists(input.seq_illumina) else "WARN", "seqkit Illumina summary presence"))
        rows.append(("seqkit_ont_present", "PASS" if os.path.exists(input.seq_ont) else "WARN", "seqkit ONT summary presence"))
        rows.append(("seqkit_assembly_present", "PASS" if os.path.exists(input.seq_asm) else "WARN", "seqkit assembly summary presence"))

        try:
            rows.append(("assembly_bp_nonzero", "PASS" if float(str(assembly_bp).replace(",", "")) > 0 else "FAIL", f"assembly_bp={assembly_bp}"))
        except Exception:
            rows.append(("assembly_bp_nonzero", "WARN", f"assembly_bp_unparseable={assembly_bp}"))

        try:
            rows.append(("assembly_contiguity", "PASS" if int(str(contigs)) >= 1 else "FAIL", f"assembly_contig_count_fasta={contigs}"))
        except Exception:
            rows.append(("assembly_contiguity", "WARN", f"assembly_contig_count_unparseable={contigs}"))

        try:
            rows.append(("n50_recorded", "PASS" if float(str(n50).replace(",", "")) > 0 else "WARN", f"N50={n50}"))
        except Exception:
            rows.append(("n50_recorded", "WARN", f"N50_unparseable={n50}"))

        try:
            rows.append(("busco_complete_recorded", "PASS" if float(busco_c) >= 0 else "WARN", f"BUSCO_complete_percent={busco_c}"))
        except Exception:
            rows.append(("busco_complete_recorded", "WARN" if RUN_BUSCO else "SKIP", f"BUSCO_complete_unparseable={busco_c}"))

        try:
            covf = float(ont_cov)
            if RUN_MODE[iso] in {"hybrid", "ont_only"}:
                rows.append(("ont_coverage", "PASS" if covf >= 20 else "WARN", f"estimated_ont_coverage_x={covf}"))
            else:
                rows.append(("ont_coverage", "NA", f"run_mode={RUN_MODE[iso]}"))
        except Exception:
            rows.append(("ont_coverage", "WARN", f"ont_coverage_unparseable={ont_cov}"))

        overall_status = "PASS"
        if any(status == "FAIL" for _, status, _ in rows):
            overall_status = "FAIL"
        elif any(status == "WARN" for _, status, _ in rows):
            overall_status = "WARN"

        rows.append(("overall_qc_status", overall_status, "aggregate QC status"))

        with open(output[0], "w", encoding="utf-8") as out:
            out.write("check\tstatus\tdetail\n")
            for check, status, detail in rows:
                out.write(f"{check}\t{status}\t{detail}\n")


rule multiqc_custom_content:
    input:
        summary=run_summary_json("{isolate}"),
        qc=qc_assessment_tsv("{isolate}"),
        versions=tool_versions_tsv("{isolate}")
    output:
        multiqc_custom_tsv("{isolate}")
    benchmark:
        bench("{isolate}", "multiqc_custom_content")
    run:
        iso = wildcards.isolate
        os.makedirs(metadata_dir(iso), exist_ok=True)

        with open(input.summary, "r", encoding="utf-8") as h:
            s = json.load(h)

        qc_map = {}
        with open(input.qc, "r", encoding="utf-8") as h:
            next(h)
            for line in h:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    qc_map[parts[0]] = {"status": parts[1], "detail": parts[2]}

        version_map = {}
        with open(input.versions, "r", encoding="utf-8") as h:
            next(h)
            for line in h:
                parts = line.rstrip("\n").split("\t", 1)
                if len(parts) == 2:
                    version_map[parts[0]] = parts[1]

        with open(output[0], "w", encoding="utf-8") as out:
            out.write(
                "sample\t"
                "run_mode\t"
                "workflow_version\t"
                "polish_mode\t"
                "executed_stages\t"
                "final_source_stage\t"
                "final_source_sha256\t"
                "target_taxid\t"
                "kraken_decision_mode\t"
                "kraken_decision\t"
                "kraken_selected_target_taxid\t"
                "kraken_selected_target_name\t"
                "kraken_selected_target_source\t"
                "kraken_decision_confidence\t"
                "kraken_awareness_state\t"
                "kraken_runner_up_name\t"
                "kraken_runner_up_percent\t"
                "kraken_selected_vs_runner_up_delta_percent\t"
                "kraken_platform_agreement_score\t"
                "assembly_bp\t"
                "assembly_contig_count_fasta\t"
                "quast_n50\t"
                "quast_gc_percent\t"
                "busco_complete_percent\t"
                "estimated_ont_coverage_x\t"
                "qc_overall\t"
                "qc_contiguity\t"
                "qc_coverage\t"
                "assembler\tinput_strategy\tlibrary_count\ttool_spades\t"
                "tool_unicycler\t"
                "tool_racon\t"
                "tool_medaka\t"
                "tool_polypolish\t"
                "tool_kraken2\t"
                "tool_busco\t"
                "tool_quast\t"
                "tool_prokka\t"
                "tool_fastqc\t"
                "tool_fastp\t"
                "tool_nanoplot\t"
                "tool_multiqc\t"
                "tool_seqkit\t"
                "tool_samtools\t"
                "tool_bwa\n"
            )

            out.write(
                f"{iso}\t"
                f"{s.get('run_mode', 'NA')}\t"
                f"{s.get('workflow_version', 'NA')}\t"
                f"{s.get('polish_mode_requested', 'NA')}\t"
                f"{','.join(s.get('executed_stages', []))}\t"
                f"{s.get('final_source_stage', 'NA')}\t"
                f"{s.get('final_source_sha256', 'NA')}\t"
                f"{s.get('target_taxid', 'NA')}\t"
                f"{s.get('kraken_decision_mode', 'NA')}\t"
                f"{s.get('kraken_decision', 'NA')}\t"
                f"{s.get('kraken_selected_target_taxid', 'NA')}\t"
                f"{s.get('kraken_selected_target_name', 'NA')}\t"
                f"{s.get('kraken_selected_target_source', 'NA')}\t"
                f"{s.get('kraken_decision_confidence', 'NA')}\t"
                f"{s.get('kraken_awareness_state', 'NA')}\t"
                f"{s.get('kraken_runner_up_name', 'NA')}\t"
                f"{s.get('kraken_runner_up_percent', 'NA')}\t"
                f"{s.get('kraken_selected_vs_runner_up_delta_percent', 'NA')}\t"
                f"{s.get('kraken_platform_agreement_score', 'NA')}\t"
                f"{s.get('assembly_bp', 'NA')}\t"
                f"{s.get('assembly_contig_count_fasta', 'NA')}\t"
                f"{s.get('quast_n50', 'NA')}\t"
                f"{s.get('quast_gc_percent', 'NA')}\t"
                f"{s.get('busco_complete_percent', 'NA')}\t"
                f"{s.get('estimated_ont_coverage_x', 'NA')}\t"
                f"{qc_map.get('overall_qc_status', {}).get('status', 'NA')}\t"
                f"{qc_map.get('assembly_contiguity', {}).get('status', 'NA')}\t"
                f"{qc_map.get('ont_coverage', {}).get('status', 'NA')}\t"
                f"{s.get('assembler', 'NA')}\t{s.get('input_strategy', 'NA')}\t{s.get('library_count', 'NA')}\t{version_map.get('spades', 'NA')}\t"
                f"{version_map.get('unicycler', 'NA')}\t"
                f"{version_map.get('racon', 'NA')}\t"
                f"{version_map.get('medaka', 'NA')}\t"
                f"{version_map.get('polypolish', 'NA')}\t"
                f"{version_map.get('kraken2', 'NA')}\t"
                f"{version_map.get('busco', 'NA')}\t"
                f"{version_map.get('quast', 'NA')}\t"
                f"{version_map.get('prokka', 'NA')}\t"
                f"{version_map.get('fastqc', 'NA')}\t"
                f"{version_map.get('fastp', 'NA')}\t"
                f"{version_map.get('nanoplot', 'NA')}\t"
                f"{version_map.get('multiqc', 'NA')}\t"
                f"{version_map.get('seqkit', 'NA')}\t"
                f"{version_map.get('samtools', 'NA')}\t"
                f"{version_map.get('bwa', 'NA')}\n"
            )

###############################################################################
# SECTION 25 — RULE ORDER
###############################################################################

ruleorder: final_copy > unicycler

# Computational records are independent of scientific outputs and survive failures.
rule computational_report:
    output:
        report(os.path.join(provenance_runtime_dir(), "computational_report.html"),
               caption="workflow/report/computational.rst", category="Computational performance")
    benchmark:
        os.path.join(provenance_runtime_dir(), "benchmarks", "computational_report.tsv")
    run:
        render_report(provenance_runtime_dir())

COMPUTATIONAL_RECORDER = install_recorder(
    workflow, BASE_DIR, os.path.join(SOURCE_ROOT, "Hybas.smk"), EFFECTIVE_CONFIG, SAMPLES,
    settings=SETTINGS)

onstart:
    if not workflow.touch:
        SETTINGS.persist()

onsuccess:
    if COMPUTATIONAL_RECORDER is not None:
        COMPUTATIONAL_RECORDER.finish("success")

onerror:
    if COMPUTATIONAL_RECORDER is not None:
        COMPUTATIONAL_RECORDER.finish("failed")

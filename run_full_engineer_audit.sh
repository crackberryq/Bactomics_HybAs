#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ============================================================
# FULL ENGINEER AUDIT WRAPPER
# Runs: S1 -> S2 -> S3 -> S4 -> S5->S6->S8/S7
# Packages everything into: ISO/reports/FINAL_AUDIT_RESULTS/
# Builds one engineer HTML report inside FINAL bundle (REPORT/index.html)
# ============================================================

ISO="${1:-isolate2}"
THREADS="${2:-8}"
PAD_CTX="${3:-2000}"   # bp context window for S6 (up/downstream)
PAD_PROOF="${4:-200}"  # bp padding for SPAdes hit extraction in S8 before S7 translation

echo "============================================================"
echo " FULL ENGINEER AUDIT"
echo " Isolate:     $ISO"
echo " Threads:     $THREADS"
echo " Context PAD: $PAD_CTX"
echo " Proof PAD:   $PAD_PROOF"
echo "============================================================"

# ---------------------------
# Paths
# ---------------------------
AUDIT_ROOT="$ISO/reports/full_audit"            # produced by run_audit.sh and S2-4 scripts
FINAL_ROOT="$ISO/reports/FINAL_AUDIT_RESULTS"   # packaged final bundle

# Ensure roots exist
mkdir -p "$AUDIT_ROOT" "$FINAL_ROOT"

# ---------------------------
# [A] Run S1 (QC steps + read stats)
# ---------------------------
./run_audit.sh "$ISO" "$THREADS"

# ---------------------------
# [B] Run S2/S3/S4 (plots & tables)
# These scripts MUST save into $ISO/reports/full_audit/*
# ---------------------------
python3 plot_results.py "$ISO" || true
python3 plot_advanced_metrics.py "$ISO" || true
python3 analyze_granularity.py "$ISO" || true

# ---------------------------
# [C] Run forensic audit (S5→S6→S8/S7) — DO NOT TOUCH LOGIC
# ---------------------------
./run_forensic_audit.sh "$ISO" "$PAD_CTX" "$PAD_PROOF"

# ---------------------------
# [D] MultiQC (audit-only)
# Restrict to audit_root/qc_steps so it doesn't scan whole isolate folder.
# ---------------------------
if command -v multiqc >/dev/null 2>&1; then
  mkdir -p "$AUDIT_ROOT/multiqc"
  multiqc "$AUDIT_ROOT/qc_steps" \
    -o "$AUDIT_ROOT/multiqc" \
    -n "${ISO}_engineer_audit_multiqc.html" \
    --force \
    >/dev/null 2>&1 || true
else
  echo "[!] multiqc not found in PATH; skipping MultiQC build."
fi

# ============================================================
# [E] PACKAGE EVERYTHING INTO ONE ENGINEER FOLDER
# ============================================================

echo
echo "============================================================"
echo " PACKAGING FINAL RESULTS -> $FINAL_ROOT"
echo "============================================================"

QC_DEST="$FINAL_ROOT/QC_ASSEMBLY/qc_steps"
FIG_DEST="$FINAL_ROOT/FIGURES"
TAB_DEST="$FINAL_ROOT/TABLES"
LOG_DEST="$FINAL_ROOT/LOGS"
MUL_DEST="$FINAL_ROOT/MULTIQC"

S5_DEST="$FINAL_ROOT/GHOST_GENES"
FOR_FLAT="$FINAL_ROOT/FORENSIC_GENES"          # flat copies (fracture/verdict txt)
FOR_FULL="$FINAL_ROOT/FORENSIC_EVIDENCE"       # full evidence folders (clickable)

REP_DEST="$FINAL_ROOT/REPORT"

mkdir -p "$QC_DEST" "$FIG_DEST" "$TAB_DEST" "$LOG_DEST" "$MUL_DEST" \
         "$S5_DEST" "$FOR_FLAT" "$FOR_FULL" "$REP_DEST"

# ---- Copy audit QC steps
if [[ -d "$AUDIT_ROOT/qc_steps" ]]; then
  rsync -a "$AUDIT_ROOT/qc_steps/" "$QC_DEST/"
fi

# ---- Copy audit tables
if [[ -d "$AUDIT_ROOT/tables" ]]; then
  rsync -a "$AUDIT_ROOT/tables/" "$TAB_DEST/" || true
fi

# ---- Copy audit figures (any ext)
if [[ -d "$AUDIT_ROOT/figures" ]]; then
  rsync -a "$AUDIT_ROOT/figures/" "$FIG_DEST/" || true
fi

# ---- Copy logs
if [[ -d "$AUDIT_ROOT/logs" ]]; then
  rsync -a "$AUDIT_ROOT/logs/" "$LOG_DEST/" || true
fi

# ---- Copy MultiQC output
if [[ -d "$AUDIT_ROOT/multiqc" ]]; then
  rsync -a "$AUDIT_ROOT/multiqc/" "$MUL_DEST/" || true
fi

# ---- Copy S5 outputs (live in ISO/reports/)
cp -f "$ISO/reports/missing_genes_recovery.csv" "$S5_DEST/" 2>/dev/null || true
cp -f "$ISO/reports/missing_targets.txt"        "$S5_DEST/" 2>/dev/null || true

# ---- Copy FULL forensic evidence folders (best for clickable browsing)
if [[ -d "$ISO/forensic_evidence" ]]; then
  rsync -a "$ISO/forensic_evidence/" "$FOR_FULL/" || true
fi

# ---- Also copy forensic outputs into a flat folder (easy scanning)
if [[ -d "$ISO/forensic_evidence" ]]; then
  for locus_dir in "$ISO/forensic_evidence"/*; do
    [[ -d "$locus_dir" ]] || continue
    locus="$(basename "$locus_dir")"

    # fracture summary(ies)
    for f in "$locus_dir"/*_fracture_summary.txt; do
      [[ -f "$f" ]] || continue
      cp -f "$f" "$FOR_FLAT/fracture_${locus}.txt"
    done

    # proof verdict
    if [[ -f "$locus_dir/proof.txt" ]]; then
      cp -f "$locus_dir/proof.txt" "$FOR_FLAT/verdict_${locus}.txt"
    fi
  done
fi

# ---------------------------
# [F] Build engineer-facing HTML report INSIDE FINAL bundle
# This is the key fix for MultiQC + evidence links.
# ---------------------------
if [[ -f "./build_engineer_report.py" ]]; then
  python3 ./build_engineer_report.py "$ISO" "$FINAL_ROOT" engineer_report.yaml || true
else
  echo "[!] build_engineer_report.py not found; skipping engineer report build."
fi

# ---------------------------
# [G] Write a small manifest.json (engineer-friendly)
# ---------------------------
python3 - <<PY
import json, os
iso="${ISO}"
final="${FINAL_ROOT}"
manifest = {
  "isolate": iso,
  "final_bundle": final,
  "how_to_open": {
    "engineer_report": f"{final}/REPORT/index.html",
    "multiqc": f"{final}/MULTIQC/{iso}_engineer_audit_multiqc.html",
    "ghost_table": f"{final}/GHOST_GENES/missing_genes_recovery.csv",
    "targets": f"{final}/GHOST_GENES/missing_targets.txt",
    "forensic_evidence_folders": f"{final}/FORENSIC_EVIDENCE/"
  },
  "parameters": {
    "threads": int("${THREADS}"),
    "PAD_CTX_bp": int("${PAD_CTX}"),
    "PAD_PROOF_bp": int("${PAD_PROOF}")
  },
  "notes": {
    "PAD_CTX_bp": "bp context extracted around the gene from the hybrid assembly (S6).",
    "PAD_PROOF_bp": "bp padding around the best SPAdes hit segment before S7 ORF integrity check (S8→S7).",
    "scope": "This bundle includes S1–S4 QC and plots plus S5–S8 forensic ghost-gene analysis."
  }
}
os.makedirs(final, exist_ok=True)
with open(os.path.join(final, "manifest.json"), "w") as f:
  json.dump(manifest, f, indent=2)
print("[+] Wrote manifest:", os.path.join(final, "manifest.json"))
PY

echo
echo "[✓] DONE"
echo "Open (single folder):"
echo "  1) $FINAL_ROOT/REPORT/index.html"
echo "  2) $FINAL_ROOT/MULTIQC/${ISO}_engineer_audit_multiqc.html"
echo "  3) $FINAL_ROOT/GHOST_GENES/missing_genes_recovery.csv"
echo "  4) $FINAL_ROOT/FORENSIC_EVIDENCE/<target>/"
echo
echo "NOTE: PAD_CTX=$PAD_CTX (context bp) and PAD_PROOF=$PAD_PROOF (SPAdes padding bp) are saved in manifest.json"

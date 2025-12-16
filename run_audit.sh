#!/bin/bash
set -euo pipefail
shopt -s nullglob

ISO=${1:-isolate2}
LINEAGE="bacillales_odb10"
THREADS=${2:-8}

AUDIT_ROOT="$ISO/reports/full_audit"
QC_DIR="$AUDIT_ROOT/qc_steps"
LOG_DIR="$AUDIT_ROOT/logs"
TAB_DIR="$AUDIT_ROOT/tables"
FIG_DIR="$AUDIT_ROOT/figures"

mkdir -p "$QC_DIR" "$LOG_DIR" "$TAB_DIR" "$FIG_DIR"

echo "========================================================"
echo " STARTING MASTER AUDIT (S1) FOR: $ISO"
echo " Output root: $AUDIT_ROOT"
echo " QC outputs:  $QC_DIR"
echo "========================================================"

run_step() {
  local name="$1"
  local file="$2"

  if [[ -f "$file" ]]; then
    echo ">>> Analyzing Assembly: $name..."
    local step_dir="$QC_DIR/$name"
    mkdir -p "$step_dir"

    # QUAST
    quast.py -o "$step_dir" "$file" --silent \
      >"$LOG_DIR/quast_${name}.log" 2>&1 || true

    # BUSCO (unique out name per step; always log)
    if ! find "$step_dir" -type f -name "short_summary*.txt" | grep -q . ; then
      local busco_run="busco_${name}"
      echo "    -> BUSCO running: $busco_run"
      busco -i "$file" -o "$busco_run" -m genome -l "$LINEAGE" \
        --out_path "$step_dir" -c "$THREADS" --force \
        >"$LOG_DIR/busco_${name}.log" 2>&1 || true
    else
      echo "    [v] BUSCO summary already exists for $name."
    fi

    echo "    [+] Done: $step_dir"
  else
    echo "    [!] Skipped: $file not found."
  fi
}

run_step "01_Unicycler_Raw" "$ISO/asm/unicycler/assembly.fasta"

for f in "$ISO"/work/assembly.racon*.fasta; do
  fname=$(basename "$f")
  round_num="${fname%.*}"
  run_step "02_Racon_${round_num##*.}" "$f"
done

run_step "03_Medaka" "$ISO/work/medaka/consensus.fasta"
run_step "04_Polypolish" "$ISO/work/assembly.polished.fasta"

echo ">>> Analyzing Read Statistics..."
STATS_FILE="$TAB_DIR/${ISO}_reads_stats.tsv"
echo -e "Stage\tFile\tFormat\tNumSeq\tSumLen\tMinLen\tAvgLen\tMaxLen" > "$STATS_FILE"

calc_stats() {
  local tag="$1"
  local pattern="$2"
  if ls $pattern >/dev/null 2>&1; then
    seqkit stats -T $pattern | awk -v t="$tag" 'NR>1{print t "\t" $0}' >> "$STATS_FILE"
  fi
}

calc_stats "Illumina_Raw"      "$ISO/raw/illumina/*.fastq.gz"
calc_stats "Illumina_Filtered" "$ISO/trimmed/*.trimmed.fq.gz"
calc_stats "ONT_Raw"           "$ISO/raw/nanopore/*.fastq.gz"
calc_stats "ONT_Filtered"      "$ISO/work/ont_filtered.fastq.gz"

echo "========================================================"
echo " AUDIT COMPLETE (S1) FOR $ISO."
echo " Outputs:"
echo "  - QC steps:   $QC_DIR/"
echo "  - Logs:       $LOG_DIR/"
echo "  - Read stats: $STATS_FILE"
echo "========================================================"

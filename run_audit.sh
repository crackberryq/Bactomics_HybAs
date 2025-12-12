cat << 'EOF' > run_audit.sh
#!/bin/bash

# 1. READ ARGUMENT (Default to isolate2 if none provided)
ISO=${1:-isolate2}
LINEAGE="bacillales_odb10"
BASE_DIR="$ISO/reports/intermediate_qc_FULL"

echo "========================================================"
echo " STARTING MASTER AUDIT FOR: $ISO"
echo " Output: $BASE_DIR"
echo "========================================================"

mkdir -p "$BASE_DIR"

# --- PART A: ASSEMBLY QC (QUAST + BUSCO) ---

run_step() {
    local name=$1
    local file=$2
    
    if [ -f "$file" ]; then
        echo ">>> Analyzing Assembly: $name..."
        mkdir -p "$BASE_DIR/$name"
        
        # Run QUAST
        quast.py -o "$BASE_DIR/$name" "$file" --silent
        
        # Run BUSCO
        # Note: We check if BUSCO was already run to save time
        if [ ! -f "$BASE_DIR/$name/short_summary.txt" ]; then
             busco -i "$file" -o "busco_out" -m genome -l "$LINEAGE" \
                  --out_path "$BASE_DIR/$name" -c 8 --force --quiet
             # Move summary to main folder for easy parsing
             cp "$BASE_DIR/$name"/busco_out/short_summary*.txt "$BASE_DIR/$name/short_summary.txt"
        else
             echo "    [v] BUSCO already exists."
        fi
        
        echo "    [+] Stats generated."
    else
        echo "    [!] Skipped: $file not found."
    fi
}

# 1. Raw Unicycler
run_step "01_Unicycler_Raw" "$ISO/asm/unicycler/assembly.fasta"

# 2. Racon (All Rounds)
for f in $ISO/work/assembly.racon*.fasta; do
    fname=$(basename "$f")
    round_num="${fname%.*}" 
    run_step "02_Racon_${round_num##*.}" "$f"
done

# 3. Medaka
run_step "03_Medaka" "$ISO/work/medaka/consensus.fasta"

# 4. Polypolish (Final)
run_step "04_Polypolish" "$ISO/work/assembly.polished.fasta"

# --- PART B: READ STATISTICS ---

echo ">>> Analyzing Read Statistics..."
STATS_FILE="$ISO/reports/reads_stats.txt"
echo "Stage,Format,NumReads,SumLen,AvgLen,Q20(%),Q30(%)" > "$STATS_FILE"

calc_stats() {
    local tag=$1
    local pattern=$2
    if ls $pattern 1> /dev/null 2>&1; then
        seqkit stats -a -T $pattern | sed '1d' | awk -v t="$tag" '{print t OFS $0}' >> "$STATS_FILE"
    fi
}

calc_stats "Illumina_Raw" "$ISO/raw/illumina/*.fastq.gz"
calc_stats "Illumina_Filtered" "$ISO/trimmed/*.trimmed.fq.gz"
calc_stats "ONT_Raw" "$ISO/raw/nanopore/*.fastq.gz"
calc_stats "ONT_Filtered" "$ISO/work/ont_filtered.fastq.gz"

echo "========================================================"
echo " AUDIT COMPLETE FOR $ISO."
echo "========================================================"
EOF
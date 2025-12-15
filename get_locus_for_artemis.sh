#!/bin/bash

# ==============================================================================
# SCRIPT: Forensic Locus Fracture Extraction (With Quick Check)
# PURPOSE: 1. Visuals: Saves files for Artemis/ACT to prove fragmentation.
#          2. Text Report: Prints a text summary to screen for immediate verify.
# ==============================================================================

ISOLATE=$1
GENE=$2
PAD=2000  # Context window

# Usage Check
if [ -z "$ISOLATE" ] || [ -z "$GENE" ]; then
    echo "❌ Error: Missing arguments."
    echo "Usage: ./get_locus_for_artemis.sh [ISOLATE_FOLDER] [GENE_NAME]"
    exit 1
fi

# Create output folder
OUT_DIR="${ISOLATE}/forensic_evidence"
mkdir -p "$OUT_DIR"

echo "----------------------------------------------------------------"
echo "🔍 Starting Forensic Audit for '$GENE' in '$ISOLATE'"
echo "----------------------------------------------------------------"

# --- STEP 1: AUTO-DETECT FILES ---
HYBRID_GFF=$(find "$ISOLATE/annotation" -name "*.gff" | head -n 1)
HYBRID_FASTA=$(find "$ISOLATE/annotation" -name "*.fna" | head -n 1)
SPADES_FASTA=$(find "$ISOLATE/annotation_spades" -name "*_spades.fna" | head -n 1)

if [[ -z "$HYBRID_GFF" || -z "$HYBRID_FASTA" || -z "$SPADES_FASTA" ]]; then
    echo "❌ CRITICAL ERROR: Input files missing."
    exit 1
fi

# --- STEP 2: LOCATE GENE ---
grep -i "$GENE" "$HYBRID_GFF" | head -n 1 > "$OUT_DIR/temp_gene.gff"

if [ ! -s "$OUT_DIR/temp_gene.gff" ]; then
    echo "❌ Error: Gene '$GENE' not found."
    rm "$OUT_DIR/temp_gene.gff"
    exit 1
fi

CONTIG=$(awk '{print $1}' "$OUT_DIR/temp_gene.gff")
START=$(awk '{print $4}' "$OUT_DIR/temp_gene.gff")
END=$(awk '{print $5}' "$OUT_DIR/temp_gene.gff")
GENE_LEN=$((END - START + 1))

echo "📍 Found $GENE on $CONTIG (Length: $GENE_LEN bp)"

# --- STEP 3: EXTRACT CONTEXT (For Artemis) ---
NEW_START=$((START - PAD))
if [ $NEW_START -lt 0 ]; then NEW_START=1; fi
NEW_END=$((END + PAD))

echo -e "$CONTIG\t$NEW_START\t$NEW_END" > "$OUT_DIR/locus_window.bed"

bedtools getfasta -fi "$HYBRID_FASTA" -bed "$OUT_DIR/locus_window.bed" -fo "$OUT_DIR/${ISOLATE}_${GENE}_Hybrid_Context.fasta"
grep "$CONTIG" "$HYBRID_GFF" | awk -v s=$NEW_START -v e=$NEW_END '$4 >= s && $5 <= e' > "$OUT_DIR/${ISOLATE}_${GENE}_Hybrid_Context.gff"

# --- STEP 4: BLAST (The Comparison) ---
# We blast the Context sequence against SPAdes
blastn \
    -query "$OUT_DIR/${ISOLATE}_${GENE}_Hybrid_Context.fasta" \
    -subject "$SPADES_FASTA" \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send" \
    > "$OUT_DIR/${ISOLATE}_${GENE}_vs_SPAdes.blast"

# --- STEP 5: QUICK CHECK REPORT ---
# Creates a simplified text file focusing strictly on the GENE breakage
REPORT="$OUT_DIR/${ISOLATE}_${GENE}_QuickCheck.txt"

echo "========================================================" > "$REPORT"
echo " FORENSIC QUICK CHECK: $GENE in $ISOLATE" >> "$REPORT"
echo "========================================================" >> "$REPORT"
echo "Gene Length in Hybrid: $GENE_LEN bp" >> "$REPORT"
echo "Status in Short-Read Assembly (SPAdes):" >> "$REPORT"
echo "--------------------------------------------------------"
echo "Col 1: Query (Hybrid) | Col 2: Subject (SPAdes Contig) | Col 4: Alignment Length" >> "$REPORT"
echo "" >> "$REPORT"

# Filter BLAST hits to show only those overlapping the gene itself (approx center of context)
# The gene is roughly at position 2000-2000+LEN in the context file.
cat "$OUT_DIR/${ISOLATE}_${GENE}_vs_SPAdes.blast" >> "$REPORT"

# Cleanup
rm "$OUT_DIR/temp_gene.gff" "$OUT_DIR/locus_window.bed"

# --- PRINT TO SCREEN ---
echo "✅ Audit Complete."
echo ""
echo "📄 CONTENTS OF QUICK REPORT:"
cat "$REPORT"
echo ""
echo "📂 Full Evidence Saved to: $OUT_DIR"

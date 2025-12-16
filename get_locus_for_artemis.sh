#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ISOLATE="${1:-}"
TARGET="${2:-}"
PAD="${3:-2000}"

if [[ -z "$ISOLATE" || -z "$TARGET" ]]; then
  echo "Usage: ./get_locus_for_artemis.sh <ISOLATE> <TARGET(gene|locus_tag|ID)> [PAD=2000]" >&2
  exit 1
fi

TARGET_CLEAN="$(printf "%s" "$TARGET" | tr -d '\r' | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"

OUT_DIR="${ISOLATE}/forensic_evidence/${TARGET_CLEAN}"
mkdir -p "$OUT_DIR"

HYBRID_GFF="$(find "${ISOLATE}/annotation" -maxdepth 1 -name "*.gff" | head -n 1)"
HYBRID_FASTA="$(find "${ISOLATE}/annotation" -maxdepth 1 -name "*.fna" -o -name "*.fa" -o -name "*.fasta" | head -n 1)"
SPADES_FASTA="${ISOLATE}/asm/spades/scaffolds.fasta"

if [[ ! -s "$HYBRID_GFF" || ! -s "$HYBRID_FASTA" || ! -s "$SPADES_FASTA" ]]; then
  echo "ERROR: missing inputs. Need:" >&2
  echo "  - $HYBRID_GFF" >&2
  echo "  - $HYBRID_FASTA" >&2
  echo "  - $SPADES_FASTA" >&2
  exit 1
fi

# --- Find the CDS feature line that matches TARGET_CLEAN ---
TARGET_LINE="$(
awk -F'\t' -v t="$TARGET_CLEAN" '
BEGIN{IGNORECASE=1}
$0 ~ /^#/ {next}
$3!="CDS" {next}
{
  # parse attributes into a map
  n=split($9, a, ";")
  id=""; locus=""; gene=""; name=""; prod=""
  for(i=1;i<=n;i++){
    gsub(/^[ \t]+|[ \t]+$/, "", a[i])
    split(a[i], kv, "=")
    k=kv[1]; v=kv[2]
    if(k=="ID") id=v
    else if(k=="locus_tag") locus=v
    else if(k=="gene") gene=v
    else if(k=="Name") name=v
    else if(k=="product") prod=v
  }

  if(id==t || locus==t || gene==t || name==t){
    print $0
    exit
  }

  # Optional: allow product substring match (comment out if you want strict matching only)
  if(prod ~ t){
    print $0
    exit
  }
}
' "$HYBRID_GFF"
)"

if [[ -z "$TARGET_LINE" ]]; then
  echo "ERROR: target '$TARGET_CLEAN' not found in $HYBRID_GFF" >&2
  exit 1
fi

CONTIG="$(printf "%s\n" "$TARGET_LINE" | awk -F'\t' '{print $1}')"
START="$( printf "%s\n" "$TARGET_LINE" | awk -F'\t' '{print $4}')"
END="$(   printf "%s\n" "$TARGET_LINE" | awk -F'\t' '{print $5}')"
STRAND="$(printf "%s\n" "$TARGET_LINE" | awk -F'\t' '{print $7}')"
LEN=$((END-START+1))

# Also extract the "Name" (gene symbol) if present, for nicer reporting
FEATURE_NAME="$(
printf "%s\n" "$TARGET_LINE" | awk -F'\t' '
{
  n=split($9,a,";"); name=""
  for(i=1;i<=n;i++){
    split(a[i],kv,"=")
    if(kv[1]=="Name"){name=kv[2]}
  }
  print name
}'
)"

if [[ -n "$FEATURE_NAME" ]]; then
  echo "Found $TARGET_CLEAN (Name=$FEATURE_NAME) on $CONTIG:${START}-${END} ($STRAND) len=${LEN} bp"
else
  echo "Found $TARGET_CLEAN on $CONTIG:${START}-${END} ($STRAND) len=${LEN} bp"
fi

# --- Extract gene-only sequence from hybrid annotation FASTA ---
echo -e "${CONTIG}\t$((START-1))\t${END}\t${TARGET_CLEAN}\t0\t${STRAND}" > "$OUT_DIR/hybrid_gene.bed"
bedtools getfasta -fi "$HYBRID_FASTA" -bed "$OUT_DIR/hybrid_gene.bed" -s -fo "$OUT_DIR/${TARGET_CLEAN}_hybrid_gene.fna"

# --- Extract context (PAD bp on both sides) from hybrid ---
WSTART=$((START-PAD)); ((WSTART<1)) && WSTART=1
WEND=$((END+PAD))
echo -e "${CONTIG}\t$((WSTART-1))\t${WEND}\t${TARGET_CLEAN}_ctx\t0\t+" > "$OUT_DIR/hybrid_ctx.bed"
bedtools getfasta -fi "$HYBRID_FASTA" -bed "$OUT_DIR/hybrid_ctx.bed" -fo "$OUT_DIR/${TARGET_CLEAN}_hybrid_context_${PAD}.fna"

# --- BLAST gene against SPAdes scaffolds to see fracture ---
blastn -query "$OUT_DIR/${TARGET_CLEAN}_hybrid_gene.fna" -subject "$SPADES_FASTA" \
  -outfmt "6 sseqid pident length qlen qstart qend sstart send bitscore" \
  | sort -k9,9nr > "$OUT_DIR/${TARGET_CLEAN}_vs_SPAdes_gene.blast"

# --- Summary report ---
REPORT="$OUT_DIR/${TARGET_CLEAN}_fracture_summary.txt"
{
  echo "Target: $TARGET_CLEAN"
  [[ -n "$FEATURE_NAME" ]] && echo "Name:   $FEATURE_NAME"
  echo "Hybrid: $CONTIG:$START-$END ($STRAND) len=$LEN bp"
  echo
  echo "Top BLAST hits (hybrid gene vs SPAdes scaffolds):"
  echo "sseqid  pident  aln_len  qlen  qstart-qend  sstart-send  bitscore"
  head -n 20 "$OUT_DIR/${TARGET_CLEAN}_vs_SPAdes_gene.blast" | awk '{print $1,$2,$3,$4,$5"-"$6,$7"-"$8,$9}'
  echo
  echo "Contigs hit (unique):"
  awk '{print $1}' "$OUT_DIR/${TARGET_CLEAN}_vs_SPAdes_gene.blast" | sort -u | nl -ba
} > "$REPORT"

echo "Saved: $REPORT"

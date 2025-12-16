#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ISOLATE="${1:-}"
TARGET="${2:-}"
PAD="${3:-200}"

if [[ -z "$ISOLATE" || -z "$TARGET" ]]; then
  echo "Usage: ./run_gene_proof.sh <ISOLATE> <TARGET(gene|locus_tag|ID)> [PAD=200]" >&2
  exit 1
fi

TARGET_CLEAN="$(printf "%s" "$TARGET" | tr -d '\r' | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"

OUT_DIR="${ISOLATE}/forensic_evidence/${TARGET_CLEAN}"
mkdir -p "$OUT_DIR"

HYBRID_GFF="$(find "${ISOLATE}/annotation" -maxdepth 1 -name "*.gff" | head -n 1)"
HYBRID_FASTA="$(find "${ISOLATE}/annotation" -maxdepth 1 -name "*.fna" -o -name "*.fa" -o -name "*.fasta" | head -n 1)"
SPADES_FASTA="${ISOLATE}/asm/spades/scaffolds.fasta"

if [[ ! -s "$HYBRID_GFF" || ! -s "$HYBRID_FASTA" || ! -s "$SPADES_FASTA" ]]; then
  echo "Missing input files." >&2
  exit 1
fi

# Find target CDS line (same robust logic as S6)
TARGET_LINE="$(
awk -F'\t' -v t="$TARGET_CLEAN" '
BEGIN{IGNORECASE=1}
$0 ~ /^#/ {next}
$3!="CDS" {next}
{
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
  if(id==t || locus==t || gene==t || name==t){ print $0; exit }
  if(prod ~ t){ print $0; exit }
}
' "$HYBRID_GFF"
)"

if [[ -z "$TARGET_LINE" ]]; then
  echo "WARN: target '$TARGET_CLEAN' not found in GFF (ID/locus_tag/gene/Name/product). Skipping." >&2
  exit 0
fi

H_CONTIG="$(printf "%s\n" "$TARGET_LINE" | awk -F'\t' '{print $1}')"
H_START="$( printf "%s\n" "$TARGET_LINE" | awk -F'\t' '{print $4}')"
H_END="$(   printf "%s\n" "$TARGET_LINE" | awk -F'\t' '{print $5}')"
H_STRAND="$(printf "%s\n" "$TARGET_LINE" | awk -F'\t' '{print $7}')"

echo "Target: $TARGET_CLEAN matched CDS on $H_CONTIG:$H_START-$H_END strand=$H_STRAND"

# Extract hybrid gene CDS (authoritative reference)
echo -e "${H_CONTIG}\t$((H_START-1))\t${H_END}\t${TARGET_CLEAN}\t0\t${H_STRAND}" > "$OUT_DIR/hybrid_gene.bed"
bedtools getfasta -fi "$HYBRID_FASTA" -bed "$OUT_DIR/hybrid_gene.bed" -s -fo "$OUT_DIR/hybrid_gene.fna"

# BLAST against SPAdes scaffolds to locate best segment
blastn -query "$OUT_DIR/hybrid_gene.fna" -subject "$SPADES_FASTA" \
  -outfmt "6 sseqid pident length qlen qstart qend sstart send bitscore" \
  | awk '$2>=90 && $3>=50 {print}' | sort -k9,9nr > "$OUT_DIR/spades_hits.tsv"

BEST="$(head -n 1 "$OUT_DIR/spades_hits.tsv" || true)"
if [[ -z "$BEST" ]]; then
  echo "No suitable BLAST hit in SPAdes for $TARGET_CLEAN (>=90% and >=50bp)."
  exit 0
fi

S_CONTIG="$(echo "$BEST" | awk '{print $1}')"
S_START="$( echo "$BEST" | awk '{print $7}')"
S_END="$(   echo "$BEST" | awk '{print $8}')"

# Normalize coordinates and strand for extraction
if (( S_START > S_END )); then
  S_STRAND="-"
  L=$S_END; R=$S_START
else
  S_STRAND="+"
  L=$S_START; R=$S_END
fi

BED_START=$((L-1-PAD)); ((BED_START<0)) && BED_START=0
BED_END=$((R+PAD))

echo -e "${S_CONTIG}\t${BED_START}\t${BED_END}\tSPAdes_hit\t0\t${S_STRAND}" > "$OUT_DIR/spades_hit.bed"
bedtools getfasta -fi "$SPADES_FASTA" -bed "$OUT_DIR/spades_hit.bed" -s -fo "$OUT_DIR/spades_locus_padded.fna"

# Run S7 proof
python3 prove_frameshift.py "$OUT_DIR/hybrid_gene.fna" "$OUT_DIR/spades_locus_padded.fna" > "$OUT_DIR/proof.txt"
echo "Saved proof: $OUT_DIR/proof.txt"

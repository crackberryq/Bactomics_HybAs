#!/bin/bash
ISOLATE=$1
GENE=$2

echo "----------------------------------------------------"
echo "🧬 PROVING GENE INTEGRITY: $GENE in $ISOLATE"
echo "----------------------------------------------------"

# 1. Extract Hybrid Gene (Reference)
HYBRID_GFF=$(find "$ISOLATE/annotation" -name "*.gff" | head -n 1)
HYBRID_FASTA=$(find "$ISOLATE/annotation" -name "*.fna" | head -n 1)

grep -i "$GENE" "$HYBRID_GFF" | head -n 1 > temp_target.gff
if [ ! -s temp_target.gff ]; then echo "❌ Gene not found"; exit 1; fi

read H_CONTIG _ _ H_START H_END _ H_STRAND _ <<< $(awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' temp_target.gff)
echo "🎯 Target: $GENE on $H_CONTIG ($H_START-$H_END) Strand: $H_STRAND"

echo -e "$H_CONTIG\t$((H_START-1))\t$H_END\t$GENE\t0\t$H_STRAND" > hybrid_gene.bed
bedtools getfasta -fi "$HYBRID_FASTA" -bed hybrid_gene.bed -s -fo hybrid_pure.fasta

# 2. Find in SPAdes
SPADES_FASTA=$(find "$ISOLATE/annotation_spades" -name "*_spades.fna" | head -n 1)
blastn -query hybrid_pure.fasta -subject "$SPADES_FASTA" -outfmt "6 sseqid sstart send length" | head -n 1 > spades_hit.txt
read S_CONTIG S_START S_END ALIGN_LEN <<< $(cat spades_hit.txt)

echo "📍 Found in SPAdes: $S_CONTIG ($S_START-$S_END)"

# 3. Extract SPAdes (Strand Logic)
if [ "$S_START" -gt "$S_END" ]; then
    BED_START=$((S_END-1)); BED_END=$S_START; S_STRAND="-"
else
    BED_START=$((S_START-1)); BED_END=$S_END; S_STRAND="+"
fi

echo -e "$S_CONTIG\t$BED_START\t$BED_END\tSPAdes_Locus\t0\t$S_STRAND" > spades_locus.bed
bedtools getfasta -fi "$SPADES_FASTA" -bed spades_locus.bed -s -fo spades_pure.fasta

# 4. Run Python Proof
python3 prove_frameshift.py hybrid_pure.fasta spades_pure.fasta

# Cleanup
rm temp_target.gff hybrid_gene.bed hybrid_pure.fasta spades_hit.txt spades_locus.bed spades_pure.fasta

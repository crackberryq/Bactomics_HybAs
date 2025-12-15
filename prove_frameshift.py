import sys
import os

# CODING TABLE (Standard Genetic Code)
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def translate_dna(dna):
    protein = ""
    if len(dna) % 3 == 0:
        for i in range(0, len(dna), 3):
            codon = dna[i:i+3].upper()
            protein += CODON_TABLE.get(codon, 'X')
    return protein

def read_fasta(filepath):
    """Reads the first sequence from a FASTA file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    seq = "".join([line.strip() for line in lines if not line.startswith(">")])
    return seq

def main():
    if len(sys.argv) < 3:
        print("Usage: python prove_frameshift.py [HYBRID_GENE_FASTA] [SPADES_LOCUS_FASTA]")
        sys.exit(1)

    hybrid_file = sys.argv[1]
    spades_file = sys.argv[2] # We assume this contains the BLAST hit region

    print(f"--- 🧬 FRAMESHIFT FORENSICS ---")
    
    # 1. Load Sequences
    hybrid_seq = read_fasta(hybrid_file)
    spades_seq = read_fasta(spades_file)
    
    print(f"Hybrid Gene Length: {len(hybrid_seq)} bp")
    print(f"SPAdes Locus Length: {len(spades_seq)} bp")

    # 2. Calculate Indel
    len_diff = len(spades_seq) - len(hybrid_seq)
    print(f"\n[MEASUREMENT]")
    print(f"Length Difference (Indel): {len_diff} bp")
    
    if len_diff % 3 != 0:
        print(f"⚠️  CRITICAL: Indel is NOT a multiple of 3. FRAMESHIFT CONFIRMED.")
    elif len_diff == 0:
        print(f"✅  Length is identical. Checking for Point Mutations...")
    else:
        print(f"ℹ️  Indel is a multiple of 3 (Codon insertion/deletion). Frame preserved.")

    # 3. Translate Comparison
    print(f"\n[BIOLOGICAL PROOF]")
    prot_hybrid = translate_dna(hybrid_seq)
    
    # Naive translation of SPAdes (assuming starts at same pos)
    # Ideally we'd align, but for a quick check, we check for early stops
    prot_spades = ""
    try:
        # Translate frame 0, 1, 2 to find best match
        best_prot = ""
        min_stops = 9999
        for frame in range(3):
            p = ""
            for i in range(frame, len(spades_seq)-2, 3):
                p += CODON_TABLE.get(spades_seq[i:i+3].upper(), 'X')
            stops = p.count('_')
            if stops < min_stops:
                min_stops = stops
                best_prot = p
        prot_spades = best_prot
    except:
        prot_spades = "Error"

    stops_hybrid = prot_hybrid.count('_')
    stops_spades = prot_spades.count('_') # stop codon is _

    print(f"Hybrid Protein Stop Codons: {stops_hybrid}")
    print(f"SPAdes Protein Stop Codons: {stops_spades}")

    if stops_spades > 1 and stops_hybrid <= 1:
        print(f"🚨 VERDICT: PREMATURE STOP CODONS DETECTED IN SPADES.")
        print(f"   The gene is broken in the short-read assembly.")
    elif stops_spades == stops_hybrid:
        print(f"✅ VERDICT: Protein appears intact. It might be a promoter/start codon mutation.")
    else:
        print(f"❓ VERDICT: Complex error. Visual inspection required.")

if __name__ == "__main__":
    main()

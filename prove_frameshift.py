#!/usr/bin/env python3
import sys

CODON_TABLE = {
   'ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T',
   'AAC':'N','AAT':'N','AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R',
   'CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P',
   'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R',
   'GTA':'V','GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A',
   'GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G',
   'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L',
   'TAC':'Y','TAT':'Y','TAA':'_','TAG':'_','TGC':'C','TGT':'C','TGA':'_','TGG':'W',
}

def read_fasta(path: str) -> str:
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq).upper()

def translate(seq: str, frame: int=0) -> str:
    prot = []
    for i in range(frame, len(seq)-2, 3):
        prot.append(CODON_TABLE.get(seq[i:i+3], 'X'))
    return "".join(prot)

def main():
    if len(sys.argv) < 3:
        print("Usage: prove_frameshift.py <HYBRID_GENE.fna> <SPADES_SEGMENT.fna>")
        sys.exit(1)

    ref = read_fasta(sys.argv[1])
    seg = read_fasta(sys.argv[2])

    print("--- FRAMESHIFT FORENSICS (S7 final) ---")
    print(f"Hybrid CDS length:   {len(ref)} bp")
    print(f"SPAdes segment len:  {len(seg)} bp")

    delta = len(seg) - len(ref)
    print(f"ΔL (seg - ref):      {delta} bp")

    mod3_ok = (delta % 3 == 0)
    print(f"Indel modulo-3:      {'OK' if mod3_ok else 'FAIL (frameshift-compatible length change)'}")

    ref_p = translate(ref, 0)
    ref_stops = ref_p.count('_')
    print(f"Ref stop codons:     {ref_stops}")

    # Best frame in segment = the one with fewest stops, then fewest X
    best = None
    for frame in (0,1,2):
        p = translate(seg, frame)
        stops = p.count('_')
        xs = p.count('X')
        score = (stops, xs)
        if best is None or score < best[0]:
            best = (score, frame, p)

    (_, _), best_frame, best_p = best
    seg_stops = best_p.count('_')
    seg_x = best_p.count('X')

    print(f"Seg best frame:      {best_frame}")
    print(f"Seg stop codons:     {seg_stops}")
    print(f"Seg X codons:        {seg_x}")

    verdict = "ORF-LIKELY-INTACT (annotation failure likely due to start/RBS/edge effects)."
    if (not mod3_ok) and (seg_stops > ref_stops):
        verdict = "FRAMESHIFT-LIKELY (non-triplet length change + excess stops)."
    elif (not mod3_ok):
        verdict = "FRAMESHIFT-POSSIBLE (non-triplet length change)."
    elif seg_stops > max(ref_stops + 1, 2):
        verdict = "BROKEN-ORF-LIKELY (excess premature stops)."

    print(f"VERDICT: {verdict}")

if __name__ == "__main__":
    main()

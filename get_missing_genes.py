import os
import sys
import pandas as pd

# --- CONFIGURATION ---
# Default to 'isolate3' if no argument provided
ISO = sys.argv[1] if len(sys.argv) > 1 else "isolate3"

# Paths to Prokka TSV files (Adjust if your folder structure differs)
HYBRID_TSV = f"{ISO}/annotation/{ISO}.tsv"
# We assume standard output path for SPAdes annotation
SHORT_READ_TSV = f"{ISO}/annotation_spades/{ISO}_spades.tsv" 
OUTPUT_CSV = f"{ISO}/reports/missing_genes_recovery.csv"

print(f"--- ANALYZING GENE RECOVERY FOR {ISO} ---")

# Check files
if not os.path.exists(HYBRID_TSV):
    print(f"ERROR: Hybrid annotation not found at {HYBRID_TSV}")
    sys.exit(1)

if not os.path.exists(SHORT_READ_TSV):
    print(f"ERROR: Short-read annotation not found at {SHORT_READ_TSV}")
    print("Run this command first: prokka --outdir {ISO}/annotation_spades --prefix {ISO}_spades {ISO}/asm/spades/scaffolds.fasta")
    sys.exit(1)

# --- 1. LOAD DATA ---
def load_genes(tsv_path):
    # Prokka TSV columns: locus_tag, ftype, length_bp, gene, EC_number, COG, product
    try:
        df = pd.read_csv(tsv_path, sep='\t')
        # Filter for CDS only
        df = df[df['ftype'] == 'CDS']
        return df
    except Exception as e:
        print(f"Error reading {tsv_path}: {e}")
        sys.exit(1)

hyb_df = load_genes(HYBRID_TSV)
short_df = load_genes(SHORT_READ_TSV)

print(f"Hybrid Total CDS:     {len(hyb_df)}")
print(f"Short-read Total CDS: {len(short_df)}")

# --- 2. FIND MISSING GENES ---
# Create reference sets from short read data
short_genes = set(short_df['gene'].dropna())
short_products = set(short_df['product'].dropna())

missing_list = []

for index, row in hyb_df.iterrows():
    gene = row['gene']
    prod = row['product']
    
    is_missing = False
    
    # Logic: 
    # 1. If gene name exists (e.g. 'ureC'), check if it's in short read set.
    # 2. If gene name is empty, check 'product' description (excluding generic "hypothetical protein").
    
    if pd.notna(gene) and gene != '':
        if gene not in short_genes:
            is_missing = True
    else:
        # Fallback to product comparison
        if pd.notna(prod) and "hypothetical protein" not in prod:
            if prod not in short_products:
                is_missing = True

    if is_missing:
        missing_list.append({
            'Locus_Tag': row['locus_tag'],
            'Gene': gene if pd.notna(gene) else "N/A",
            'Product': prod,
            'Length_bp': row['length_bp']
        })

# --- 3. SAVE RESULTS ---
if missing_list:
    miss_df = pd.DataFrame(missing_list)
    print(f"\n[!] FOUND {len(miss_df)} GENES UNIQUE TO HYBRID ASSEMBLY!")
    print(miss_df.head(10)) # Show preview
    miss_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\n[+] Saved list to: {OUTPUT_CSV}")
else:
    print("\n[?] No unique non-hypothetical genes found.")

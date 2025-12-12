cat << 'EOF' > analyze_granularity.py
import os
import sys
import glob
import re
import subprocess
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# --- 1. GET ARGUMENT ---
if len(sys.argv) < 2:
    print("Usage: python analyze_granularity.py <isolate_name>")
    print("Example: python analyze_granularity.py isolate2")
    sys.exit(1)

ISO = sys.argv[1]
QC_DIR = f"{ISO}/reports/intermediate_qc_FULL"
OUTPUT_DIR = f"{ISO}/reports"

# Map step names to actual files for dnadiff
FILE_MAP = {
    "01_Unicycler_Raw": f"{ISO}/asm/unicycler/assembly.fasta",
    "02_Racon_racon0":  f"{ISO}/work/assembly.racon0.fasta", 
    "02_Racon_racon1":  f"{ISO}/work/assembly.racon1.fasta",
    "02_Racon_racon2":  f"{ISO}/work/assembly.racon2.fasta",
    "03_Medaka":        f"{ISO}/work/medaka/consensus.fasta",
    "04_Polypolish":    f"{ISO}/work/assembly.polished.fasta"
}

print(f"--- ANALYZING {ISO} ---")

if not os.path.exists(QC_DIR):
    print(f"ERROR: {QC_DIR} not found. Run bash run_audit.sh {ISO} first.")
    sys.exit(1)

data = []
steps = sorted([d for d in os.listdir(QC_DIR) if os.path.isdir(os.path.join(QC_DIR, d))])

prev_fasta = None

for step in steps:
    path = os.path.join(QC_DIR, step)
    
    # A. Get N50
    n50 = 0
    if os.path.exists(f"{path}/report.txt"):
        with open(f"{path}/report.txt") as f:
            for line in f:
                if line.startswith("N50"):
                    try: n50 = int(line.split()[1])
                    except: pass
                    break

    # B. Get BUSCO
    comp = 0.0
    frag = 0.0
    b_files = glob.glob(f"{path}/**/short_summary*.txt", recursive=True)
    if b_files:
        with open(b_files[0]) as f:
            for line in f:
                if "C:" in line:
                    try:
                        comp = float(re.search(r'C:\s*([\d\.]+)%', line).group(1))
                        frag = float(re.search(r'F:\s*([\d\.]+)%', line).group(1))
                    except: pass
                    break

    # C. Get Deltas (dnadiff)
    snps = 0
    indels = 0
    current_fasta = None
    
    # Find matching file for this step
    for key in FILE_MAP:
        if key in step: 
            current_fasta = FILE_MAP[key]
            break
            
    if prev_fasta and current_fasta and os.path.exists(prev_fasta) and os.path.exists(current_fasta):
        subprocess.run(f"dnadiff -p temp {prev_fasta} {current_fasta}", shell=True, stderr=subprocess.DEVNULL)
        if os.path.exists("temp.report"):
            with open("temp.report") as f:
                for line in f:
                    if "TotalSNPs" in line: snps = int(line.split()[1])
                    if "TotalIndels" in line: indels = int(line.split()[1])
            os.system("rm temp.*")

    if current_fasta and os.path.exists(current_fasta):
        prev_fasta = current_fasta

    # Clean Label
    label = step.split('_', 1)[1].replace('_', ' ')
    label = label.replace("Racon racon", "Racon Round")
    
    data.append({'Step': label, 'N50': n50, 'Complete': comp, 'Fragmented': frag, 'SNPs': snps, 'Indels': indels})

df = pd.DataFrame(data)
print(df)
df.to_csv(f"{OUTPUT_DIR}/{ISO}_granular_stats.csv", index=False)

# PLOTTING
sns.set_theme(style="whitegrid")
fig = plt.figure(figsize=(14, 8))
gs = fig.add_gridspec(2, 1, height_ratios=[2, 1])

# Plot 1: BUSCO Recovery
ax1 = fig.add_subplot(gs[0])
ax1.plot(df['Step'], df['Complete'], marker='o', linewidth=3, color='#2ca02c', label='BUSCO Complete %')
ax1.set_ylabel('Completeness (%)', fontweight='bold', color='#2ca02c')
ax1.set_title(f'Pipeline Evolution: {ISO}', fontweight='bold', fontsize=14)
ax1.legend(loc='lower right')
ax1.grid(True, alpha=0.3)

# Annotate N50 on the same plot
ax1b = ax1.twinx()
ax1b.plot(df['Step'], df['N50'], marker='s', linestyle='--', color='#1f77b4', label='N50 (bp)', alpha=0.6)
ax1b.set_ylabel('N50 (bp)', fontweight='bold', color='#1f77b4')

# Plot 2: Corrections
ax2 = fig.add_subplot(gs[1])
x = range(len(df))
ax2.bar(x, df['Indels'], color='#e67e22', label='Indels Fixed')
ax2.bar(x, df['SNPs'], bottom=df['Indels'], color='#e74c3c', label='SNPs Fixed')
ax2.set_xticks(x)
ax2.set_xticklabels(df['Step'], rotation=45, ha='right')
ax2.set_ylabel('Modifications vs Prev', fontweight='bold')
ax2.legend()

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/{ISO}_evolution.png", dpi=300)
print(f"[+] Saved: {OUTPUT_DIR}/{ISO}_evolution.png")
EOF
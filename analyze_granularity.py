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
    sys.exit(1)

ISO = sys.argv[1]
QC_DIR = f"{ISO}/reports/intermediate_qc_FULL"
OUTPUT_DIR = f"{ISO}/reports"

# Map step names to actual files for dnadiff
# Adjust paths if your specific file naming differs
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
df.to_csv(f"{OUTPUT_DIR}/{ISO}_granular_stats.csv", index=False)

# --- PLOTTING (THE 3-PANEL DASHBOARD) ---
sns.set_theme(style="whitegrid")
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 2)

# Panel 1: The "U-Shape" Recovery (BUSCO) - Top Row Full Width
ax1 = fig.add_subplot(gs[0, :])
color_busco = '#2ca02c'
ax1.plot(df['Step'], df['Complete'], marker='o', linewidth=3, color=color_busco, label='BUSCO Complete %')
ax1.set_ylabel('BUSCO Completeness (%)', color=color_busco, fontweight='bold', fontsize=12)
ax1.tick_params(axis='y', labelcolor=color_busco)
ax1.set_title(f'Biological Recovery Trajectory ({ISO})', fontweight='bold', fontsize=14)
ax1.grid(True, alpha=0.3)

# Add Narrative Annotation
if not df.empty:
    min_val = df['Complete'].min()
    min_idx = df['Complete'].idxmin()
    # Annotate the lowest point
    ax1.annotate('Indel Damage\n(Long-read polishing)', 
                 xy=(min_idx, min_val), xytext=(min_idx, min_val+3),
                 arrowprops=dict(facecolor='black', shrink=0.05), ha='center', fontsize=11)

# Panel 2: Structural Stability (N50) - Bottom Left
ax2 = fig.add_subplot(gs[1, 0])
color_n50 = '#1f77b4'
ax2.plot(df['Step'], df['N50'], marker='s', linestyle='--', color=color_n50)
ax2.set_ylabel('N50 (bp)', fontweight='bold')
ax2.set_title('Structural Stability (N50)', fontweight='bold')
ax2.tick_params(axis='x', rotation=45)

# Panel 3: Granular Corrections - Bottom Right
ax3 = fig.add_subplot(gs[1, 1])
x = range(len(df))
ax3.bar(x, df['Indels'], color='#e67e22', alpha=0.8, label='Indels Modified')
ax3.bar(x, df['SNPs'], bottom=df['Indels'], color='#e74c3c', alpha=0.8, label='SNPs Modified')
ax3.set_xticks(x)
ax3.set_xticklabels(df['Step'], rotation=45, ha='right')
ax3.set_ylabel('Count (vs Previous Step)', fontweight='bold')
ax3.set_title('Granular Corrections per Step', fontweight='bold')
ax3.legend()

plt.tight_layout()
final_plot = f"{OUTPUT_DIR}/{ISO}_granular_evolution.jpg" # Saving as JPG to match your preference
plt.savefig(final_plot, dpi=300)
print(f"\n[+] Visualization saved to: {final_plot}")

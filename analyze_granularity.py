#!/usr/bin/env python3
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

# --- 1. SETUP ---
if len(sys.argv) < 2:
    print("Usage: python analyze_granularity.py <isolate_name>")
    sys.exit(1)

ISO = sys.argv[1]
AUDIT_ROOT = f"{ISO}/reports/full_audit"
QC_DIR     = f"{AUDIT_ROOT}/qc_steps"
TABLE_DIR  = f"{AUDIT_ROOT}/tables"
FIG_DIR    = f"{AUDIT_ROOT}/figures"

os.makedirs(TABLE_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# Map step names to actual files for dnadiff
FILE_MAP = {
    "01_Unicycler_Raw": f"{ISO}/asm/unicycler/assembly.fasta",
    "02_Racon_racon0":  f"{ISO}/work/assembly.racon0.fasta",
    "02_Racon_racon1":  f"{ISO}/work/assembly.racon1.fasta",
    "02_Racon_racon2":  f"{ISO}/work/assembly.racon2.fasta",
    "03_Medaka":        f"{ISO}/work/medaka/consensus.fasta",
    "04_Polypolish":    f"{ISO}/work/assembly.polished.fasta"
}

print(f"--- ANALYZING GRANULARITY FOR {ISO} ---")

if not os.path.isdir(QC_DIR):
    print(f"ERROR: QC_DIR not found: {QC_DIR}\nRun run_audit.sh first.")
    sys.exit(1)

# --- 2. DATA EXTRACTION ---
data = []
steps = sorted([d for d in os.listdir(QC_DIR) if os.path.isdir(os.path.join(QC_DIR, d))])
prev_fasta = None

for step in steps:
    path = os.path.join(QC_DIR, step)

    # A. Get N50 (Structure)
    n50 = 0
    rpt = os.path.join(path, "report.txt")
    if os.path.exists(rpt):
        with open(rpt) as f:
            for line in f:
                if line.startswith("N50"):
                    try: n50 = int(line.split()[1])
                    except: pass
                    break

    # B. Get BUSCO (Biology)
    comp, frag = 0.0, 0.0
    b_files = glob.glob(f"{path}/**/short_summary*.txt", recursive=True)
    if b_files:
        with open(b_files[0]) as f:
            content = f.read()
            mC = re.search(r'C:\s*([\d\.]+)%', content)
            mF = re.search(r'F:\s*([\d\.]+)%', content)
            if mC: comp = float(mC.group(1))
            if mF: frag = float(mF.group(1))

    # C. Get Deltas (dnadiff)
    snps, indels = 0, 0
    current_fasta = None
    for key in FILE_MAP:
        if key in step:
            current_fasta = FILE_MAP[key]
            break

    if prev_fasta and current_fasta and os.path.exists(prev_fasta) and os.path.exists(current_fasta):
        prefix = f"temp_diff_{step}"
        subprocess.run(f"dnadiff -p {prefix} {prev_fasta} {current_fasta}", 
                       shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        if os.path.exists(f"{prefix}.report"):
            with open(f"{prefix}.report") as f:
                for line in f:
                    if "TotalSNPs" in line:
                        try: snps = int(line.split()[1])
                        except: pass
                    if "TotalIndels" in line:
                        try: indels = int(line.split()[1])
                        except: pass
            subprocess.run(f"rm -f {prefix}.*", shell=True)

    if current_fasta and os.path.exists(current_fasta):
        prev_fasta = current_fasta

    # Clean Label
    label = step.split('_', 1)[1].replace('_', ' ')
    label = label.replace("Racon racon", "Racon Rd")
    
    data.append({
        'Step': label,
        'N50': n50,
        'Complete': comp,
        'Fragmented': frag,
        'SNPs': snps,
        'Indels': indels
    })

df = pd.DataFrame(data)
out_csv = f"{TABLE_DIR}/{ISO}_granular_stats.csv"
df.to_csv(out_csv, index=False)
print(f"[+] Table saved: {out_csv}")

# --- 3. PROFESSIONAL PLOTTING ---
if df.empty:
    print("[!] No data to plot.")
    sys.exit(0)

sns.set_theme(style="whitegrid")
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 2)

# PANEL 1: Biological Recovery (BUSCO) - Top Row
ax1 = fig.add_subplot(gs[0, :])
color_bio = '#2ca02c' # Nature Green
ax1.plot(df['Step'], df['Complete'], marker='o', markersize=8, linewidth=3, 
         color=color_bio, label='BUSCO Complete %')

ax1.set_ylabel('BUSCO Completeness (%)', fontweight='bold', color=color_bio)
ax1.set_title(f'Biological Recovery Trajectory ({ISO})', fontweight='bold', fontsize=14, pad=15)
ax1.grid(True, alpha=0.3)
ax1.tick_params(axis='y', labelcolor=color_bio)

# Annotation: Highlight the "Dip"
if not df.empty:
    min_idx = df['Complete'].idxmin()
    min_val = df.loc[min_idx, 'Complete']
    if min_val < df['Complete'].iloc[0] and min_val < df['Complete'].iloc[-1]:
        ax1.annotate('Indel Damage\n(Long-read polishing)', 
                     xy=(min_idx, min_val), 
                     xytext=(min_idx, min_val + 5), 
                     arrowprops=dict(facecolor='black', shrink=0.05),
                     ha='center', fontsize=11, fontweight='bold', color='#c0392b')

# PANEL 2: Structural Stability (N50) - Bottom Left
ax2 = fig.add_subplot(gs[1, 0])
color_str = '#1f77b4' # Scientific Blue
ax2.plot(df['Step'], df['N50'], marker='s', linestyle='--', linewidth=2, color=color_str)

ax2.set_ylabel('N50 (bp)', fontweight='bold', color=color_str)
ax2.set_title('Structural Stability', fontweight='bold', fontsize=12)
ax2.tick_params(axis='x', rotation=30)
ax2.grid(True, alpha=0.3)

# REVERTED: Standard notation (with 1e6 offset) is better for small ranges
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# PANEL 3: Granular Corrections (Stacked Bar) - Bottom Right
ax3 = fig.add_subplot(gs[1, 1])
x = range(len(df))
p1 = ax3.bar(x, df['Indels'], color='#e67e22', alpha=0.85, label='Indels Fixed')
p2 = ax3.bar(x, df['SNPs'], bottom=df['Indels'], color='#e74c3c', alpha=0.85, label='SNPs Fixed')

ax3.set_xticks(list(x))
ax3.set_xticklabels(df['Step'], rotation=30, ha='right')
ax3.set_ylabel('Corrections Count', fontweight='bold')
ax3.set_title('Polishing Activity per Step', fontweight='bold', fontsize=12)
ax3.legend(loc='upper right', frameon=True)
ax3.grid(axis='y', alpha=0.3)

# Final Layout
plt.tight_layout()
final_plot = f"{FIG_DIR}/{ISO}_granular_evolution.png"
plt.savefig(final_plot, dpi=300)
print(f"[+] Visualization saved: {final_plot}")
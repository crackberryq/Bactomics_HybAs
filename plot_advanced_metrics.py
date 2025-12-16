#!/usr/bin/env python3
import os
import sys
import subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# --- 1. SETUP & PATHS ---
ISO = sys.argv[1] if len(sys.argv) > 1 else "isolate2"
AUDIT_ROOT = f"{ISO}/reports/full_audit"
# Ensure the paths match your structure exactly
WORK_DIR   = f"{ISO}/work"
ANN_DIR    = f"{ISO}/annotation"
OUT_FIG    = f"{AUDIT_ROOT}/figures/{ISO}_advanced_metrics.png"

# Create output folder
os.makedirs(os.path.dirname(OUT_FIG), exist_ok=True)

# Input Files
MEDAKA_ASM = f"{WORK_DIR}/medaka/consensus.fasta"
FINAL_ASM  = f"{WORK_DIR}/assembly.polished.fasta"
BAM_FILE   = f"{WORK_DIR}/illumina.bam"
PROKKA_TXT = f"{ANN_DIR}/{ISO}.txt"

print(f"--- GENERATING ADVANCED DASHBOARD FOR {ISO} ---")

# --- 2. DATA EXTRACTION ---

# A. Polishing Corrections (dnadiff)
snps, indels = 0, 0
if os.path.exists(MEDAKA_ASM) and os.path.exists(FINAL_ASM):
    # Quietly run dnadiff
    subprocess.run(f"dnadiff -p polish_delta {MEDAKA_ASM} {FINAL_ASM}", 
                   shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if os.path.exists("polish_delta.report"):
        with open("polish_delta.report") as f:
            for line in f:
                if "TotalSNPs" in line:
                    try: snps = int(line.split()[1])
                    except: pass
                if "TotalIndels" in line:
                    try: indels = int(line.split()[1])
                    except: pass
        subprocess.run("rm -f polish_delta.*", shell=True)
else:
    print("[!] Warning: Assembly files missing for dnadiff.")

# B. Mapping Rate (Samtools)
mapped_pct = 0.0
if os.path.exists(BAM_FILE):
    # Run flagstat
    r = subprocess.run(f"samtools flagstat {BAM_FILE}", shell=True, capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if "mapped (" in line:
            try:
                # Format: 1234 + 0 mapped (99.85% : N/A)
                mapped_pct = float(line.split("(")[1].split("%")[0])
            except: pass
            break
else:
    print("[!] Warning: BAM file missing.")

# C. Biological Features (Prokka)
cds, rrna, trna = 0, 0, 0
if os.path.exists(PROKKA_TXT):
    with open(PROKKA_TXT) as f:
        for line in f:
            if "CDS:" in line:
                try: cds = int(line.split()[1])
                except: pass
            if "rRNA:" in line:
                try: rrna = int(line.split()[1])
                except: pass
            if "tRNA:" in line:
                try: trna = int(line.split()[1])
                except: pass
else:
    print("[!] Warning: Prokka summary missing.")

print(f"Stats -> SNPs: {snps}, Indels: {indels}, Mapped: {mapped_pct}%, CDS: {cds}")

# --- 3. PROFESSIONAL PLOTTING ---
sns.set_theme(style="white") # Clean canvas
fig = plt.figure(figsize=(15, 5))
gs = fig.add_gridspec(1, 3)

# PANEL A: Corrections (Vertical Bar)
ax1 = fig.add_subplot(gs[0, 0])
x_labels = ['SNPs Fixed', 'Indels Fixed']
y_values = [snps, indels]
colors_a = ['#e74c3c', '#e67e22'] # Red, Orange
bars = ax1.bar(x_labels, y_values, color=colors_a, alpha=0.85, width=0.6)

ax1.set_title(f'Polypolish Corrections', fontweight='bold', fontsize=12, pad=10)
ax1.set_ylabel('Count', fontsize=10)
ax1.bar_label(bars, padding=3, fontweight='bold', fontsize=11)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_color('#cccccc')
ax1.spines['bottom'].set_color('#cccccc')

# PANEL B: Mapping Rate (Donut Chart)
ax2 = fig.add_subplot(gs[0, 1])
if mapped_pct > 0:
    sizes = [mapped_pct, 100 - mapped_pct]
    labels = ['', ''] # No labels on outer ring
    colors_b = ['#2ecc71', '#ecf0f1'] # Green, Light Gray
    
    # Donut logic
    wedges, texts = ax2.pie(sizes, labels=labels, colors=colors_b, startangle=90, 
                            wedgeprops=dict(width=0.4, edgecolor='w'))
    
    # Center Text
    ax2.text(0, 0, f"{mapped_pct}%", ha='center', va='center', fontsize=20, fontweight='bold', color='#2c3e50')
    ax2.text(0, -0.25, "Mapped Reads", ha='center', va='center', fontsize=10, color='#7f8c8d')
    ax2.set_title('Read Consistency', fontweight='bold', fontsize=12, pad=10)
else:
    ax2.text(0.5, 0.5, "No Data", ha='center', va='center')

# PANEL C: Biology (Horizontal Bar)
ax3 = fig.add_subplot(gs[0, 2])
labels_c = ['rRNA', 'tRNA', 'CDS (x1k)']
# Normalize CDS to thousands for scale compatibility
values_c = [rrna, trna, cds/1000] 
colors_c = ['#3498db', '#9b59b6', '#34495e'] # Blue, Purple, Dark Blue

hbars = ax3.barh(labels_c, values_c, color=colors_c, alpha=0.85, height=0.6)
ax3.set_title('Biological Integrity', fontweight='bold', fontsize=12, pad=10)

# Custom labels
for i, v in enumerate(values_c):
    # formatting: if it's the CDS bar (index 2), show real count in parens
    label_text = str(int(v)) if i < 2 else f"{v:.1f} ({cds})"
    ax3.text(v + (max(values_c)*0.05), i, label_text, va='center', fontweight='bold', fontsize=10)

ax3.set_xlim(0, max(values_c)*1.3) # room for labels
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_color('#cccccc')
ax3.spines['left'].set_visible(False) # Clean look
ax3.tick_params(axis='y', length=0) # Hide ticks

# Final Layout
plt.suptitle(f"Advanced Quality Audit: {ISO}", fontsize=16, fontweight='bold', y=1.05)
plt.tight_layout()
plt.savefig(OUT_FIG, dpi=300, bbox_inches='tight')
print(f"[+] Dashboard saved: {OUT_FIG}")
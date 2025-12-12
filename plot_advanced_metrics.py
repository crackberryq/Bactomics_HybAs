import os
import sys
import subprocess
import matplotlib
matplotlib.use('Agg') # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns

# --- 1. GET ARGUMENT ---
if len(sys.argv) < 2:
    # Default to isolate3 if no name provided
    ISO = "isolate3"
else:
    ISO = sys.argv[1]

WORK_DIR = f"{ISO}/work"
ANN_DIR = f"{ISO}/annotation"
OUT_FIG = f"{ISO}/reports/{ISO}_advanced_metrics.png"

# Paths
MEDAKA_ASM = f"{WORK_DIR}/medaka/consensus.fasta"
FINAL_ASM = f"{WORK_DIR}/assembly.polished.fasta"
BAM_FILE = f"{WORK_DIR}/illumina.bam"
PROKKA_TXT = f"{ANN_DIR}/{ISO}.txt"

print(f"--- GENERATING DASHBOARD (FIG S2) FOR {ISO} ---")

# --- 1. CALCULATE POLISHING DELTAS (MUMmer) ---
snps = 0
indels = 0
if os.path.exists(MEDAKA_ASM) and os.path.exists(FINAL_ASM):
    # Quietly run dnadiff
    subprocess.run(f"dnadiff -p polish_delta {MEDAKA_ASM} {FINAL_ASM}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    if os.path.exists("polish_delta.report"):
        with open("polish_delta.report") as f:
            for line in f:
                if "TotalSNPs" in line:
                    try: snps = int(line.split()[1])
                    except: pass
                if "TotalIndels" in line:
                    try: indels = int(line.split()[1])
                    except: pass
        os.system("rm polish_delta.*")
else:
    print("[!] Warning: Assembly files missing for dnadiff.")

# --- 2. CALCULATE MAPPING RATE (Samtools) ---
mapped_pct = 0.0
if os.path.exists(BAM_FILE):
    # Run flagstat
    result = subprocess.run(f"samtools flagstat {BAM_FILE}", shell=True, capture_output=True, text=True)
    for line in result.stdout.splitlines():
        if "mapped (" in line:
            try:
                # Format: 1234 + 0 mapped (99.85% : N/A)
                mapped_pct = float(line.split('(')[1].split('%')[0])
            except: pass
            break
else:
    print("[!] Warning: BAM file missing.")

# --- 3. GET ANNOTATION COUNTS (Prokka) ---
cds = 0
rrna = 0
trna = 0
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

print(f"  Fixed:  {snps} SNPs, {indels} Indels")
print(f"  Mapped: {mapped_pct}%")

# --- 4. PLOTTING ---
sns.set_theme(style="white")
fig = plt.figure(figsize=(14, 5))
gs = fig.add_gridspec(1, 3)

# PANEL A: Corrections (Bar Chart)
ax1 = fig.add_subplot(gs[0, 0])
x_labels = ['SNPs Fixed', 'Indels Fixed']
y_values = [snps, indels]
colors = ['#e74c3c', '#e67e22'] 
bars = ax1.bar(x_labels, y_values, color=colors, alpha=0.8)
ax1.set_title(f'Polypolish Corrections\n({ISO})', fontweight='bold')
ax1.set_ylabel('Count')
ax1.bar_label(bars, padding=3, fontweight='bold')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# PANEL B: Mapping (Donut Chart)
if mapped_pct > 0:
    ax2 = fig.add_subplot(gs[0, 1])
    sizes = [mapped_pct, 100 - mapped_pct]
    labels = [f'Mapped\n({mapped_pct}%)', 'Unmapped']
    colors_pie = ['#2ecc71', '#ecf0f1']
    ax2.pie(sizes, labels=labels, colors=colors_pie, startangle=90, pctdistance=0.85, 
            explode=(0, 0.1), textprops={'fontweight': 'bold'})
    centre_circle = plt.Circle((0,0),0.70,fc='white')
    ax2.add_artist(centre_circle)
    ax2.set_title('Read Consistency\n(Illumina Mapping Rate)', fontweight='bold')

# PANEL C: Biology (Horizontal Bar)
ax3 = fig.add_subplot(gs[0, 2])
labels_bio = ['rRNA Operons', 'tRNAs', 'CDS (x1000)']
values_bio = [rrna, trna, cds/1000] 
colors_bio = ['#3498db', '#9b59b6', '#34495e']
hbars = ax3.barh(labels_bio, values_bio, color=colors_bio, alpha=0.8)
ax3.set_title('Biological Integrity\n(Features Recovered)', fontweight='bold')
if max(values_bio) > 0:
    ax3.set_xlim(0, max(values_bio)*1.2)
ax3.bar_label(hbars, fmt='%.1f', padding=3, fontweight='bold')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.text(max(values_bio)/2, 2, f"  ({cds} total)", va='center', fontsize=9, color='#34495e')

plt.tight_layout()
plt.savefig(OUT_FIG, dpi=300)
print(f"[+] Dashboard saved to: {OUT_FIG}")

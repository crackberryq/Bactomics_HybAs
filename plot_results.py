import os
import glob
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# --- CONFIGURATION ---
ISO = "isolate3"
QC_DIR = f"{ISO}/reports/intermediate_qc_FULL"
OUTPUT_FILE = f"{ISO}_benchmark_figure.png"

print(f"--- DEBUGGING DIRECTORY: {QC_DIR} ---")

data = []

# 1. Get folders
if not os.path.exists(QC_DIR):
    print("ERROR: Directory does not exist!")
    exit()
    
steps = sorted([d for d in os.listdir(QC_DIR) if os.path.isdir(os.path.join(QC_DIR, d))])

for step in steps:
    path = os.path.join(QC_DIR, step)
    
    # --- GET N50 ---
    n50 = 0
    q_file = os.path.join(path, "report.txt")
    if os.path.exists(q_file):
        with open(q_file) as f:
            for line in f:
                if line.startswith("N50"):
                    try:
                        n50 = int(line.split()[1])
                    except: pass
                    break
    
    # --- GET BUSCO (ROBUST REGEX) ---
    comp = 0.0
    frag = 0.0
    
    # Find ANY text file starting with short_summary inside the folder
    found_files = glob.glob(f"{path}/**/short_summary*.txt", recursive=True)
    
    if found_files:
        target_file = found_files[0]
        with open(target_file) as f:
            content = f.read()
            # Regex to match: "C:99.8%" OR "C: 99.8%" (handles spaces)
            match = re.search(r'C:\s*(\d+\.\d+)%.*F:\s*(\d+\.\d+)%', content)
            if match:
                comp = float(match.group(1))
                frag = float(match.group(2))
                print(f"OK: {step} -> BUSCO {comp}% (File: {os.path.basename(target_file)})")
            else:
                print(f"FAIL: {step} -> Regex did not match inside {os.path.basename(target_file)}")
                # Print a snippet to help debug
                print(f"    Snippet: {content[:100]}...") 
    else:
        print(f"MISSING: No BUSCO summary file found in {step}")

    # Label Cleaning
    label = step.split('_', 1)[1].replace('_', ' ')
    label = label.replace("Racon racon0", "Racon (Seed)")
    label = label.replace("Racon racon1", "Racon (Round 1)")
    label = label.replace("Racon racon2", "Racon (Round 2)")
    label = label.replace("Polypolish Final", "Polypolish")
    
    data.append({'Step': label, 'N50': n50, 'Complete': comp, 'Fragmented': frag})

df = pd.DataFrame(data)
print("\n--- FINAL DATA ---")
print(df)

# --- PLOT ---
if df['Complete'].max() > 0:
    sns.set_theme(style="whitegrid")
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # N50 Line
    color_n50 = '#1f77b4'
    ax1.plot(df['Step'], df['N50'], color=color_n50, marker='o', linewidth=3, label='N50 (bp)')
    ax1.set_ylabel('Assembly N50 (bp)', color=color_n50, fontweight='bold', fontsize=12)
    ax1.tick_params(axis='y', labelcolor=color_n50)
    ax1.grid(False)

    # BUSCO Bars
    ax2 = ax1.twinx()
    color_busco = '#2ca02c'
    ax2.bar(df['Step'], df['Complete'], color=color_busco, alpha=0.3, label='BUSCO %', width=0.5)
    ax2.set_ylabel('BUSCO Completeness (%)', color=color_busco, fontweight='bold', fontsize=12)
    ax2.tick_params(axis='y', labelcolor=color_busco)
    
    # Dynamic Zoom
    min_b = max(0, df['Complete'].min() - 1.0)
    ax2.set_ylim(min_b, 100.1)

    plt.title(f"Optimization Trajectory: {ISO}", fontweight='bold', fontsize=14)
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE, dpi=300)
    print(f"\n[+] Plot generated: {OUTPUT_FILE}")
else:
    print("\n[!] Skipping plot: No valid BUSCO data found.")

#!/usr/bin/env python3
import os
import glob
import re
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# --- PARSING FUNCTIONS ---
def parse_quast_n50(report_txt: str) -> int:
    if not os.path.exists(report_txt):
        return 0
    with open(report_txt) as f:
        for line in f:
            if line.startswith("N50"):
                try:
                    return int(line.split()[1])
                except:
                    return 0
    return 0

def parse_busco_summary(summary_path: str):
    # Returns: comp, frag, miss, dup, n_busco
    comp = frag = miss = dup = 0.0
    n_busco = 0
    if not summary_path or not os.path.exists(summary_path):
        return comp, frag, miss, dup, n_busco

    text = open(summary_path).read()
    
    # Modern BUSCO format
    m = re.search(r"C:\s*([\d\.]+)%.*?S:\s*([\d\.]+)%.*?D:\s*([\d\.]+)%.*?F:\s*([\d\.]+)%.*?M:\s*([\d\.]+)%.*?n:\s*(\d+)", text, re.S)
    if m:
        return float(m.group(1)), float(m.group(4)), float(m.group(5)), float(m.group(3)), int(m.group(6))

    # Fallback: older formats
    m2 = re.search(r"C:\s*([\d\.]+)%.+F:\s*([\d\.]+)%.+M:\s*([\d\.]+)%", text, re.S)
    if m2:
        comp = float(m2.group(1))
        frag = float(m2.group(2))
        miss = float(m2.group(3))
    m3 = re.search(r"n:\s*(\d+)", text)
    if m3:
        n_busco = int(m3.group(1))
    return comp, frag, miss, dup, n_busco

def find_busco_summary(step_dir: str):
    hits = glob.glob(os.path.join(step_dir, "**", "short_summary*.txt"), recursive=True)
    return hits[0] if hits else None

# --- MAIN ---
def main():
    ISO = sys.argv[1] if len(sys.argv) > 1 else "isolate2"

    # Your specific file structure
    AUDIT_ROOT = f"{ISO}/reports/full_audit"
    QC_DIR     = f"{AUDIT_ROOT}/qc_steps"
    OUT_TABLE  = f"{AUDIT_ROOT}/tables/{ISO}_benchmark_table.tsv"
    OUT_FIG    = f"{AUDIT_ROOT}/figures/{ISO}_benchmark_figure.png"

    os.makedirs(os.path.dirname(OUT_TABLE), exist_ok=True)
    os.makedirs(os.path.dirname(OUT_FIG),   exist_ok=True)

    print(f"--- BENCHMARKING {ISO} ---")
    print(f"Reading metrics from: {QC_DIR}")

    if not os.path.isdir(QC_DIR):
        print(f"ERROR: QC_DIR missing: {QC_DIR}")
        sys.exit(1)

    steps = sorted([d for d in os.listdir(QC_DIR) if os.path.isdir(os.path.join(QC_DIR, d))])

    rows = []
    for step in steps:
        step_dir = os.path.join(QC_DIR, step)

        # 1. QUAST
        n50 = parse_quast_n50(os.path.join(step_dir, "report.txt"))

        # 2. BUSCO
        summary = find_busco_summary(step_dir)
        comp, frag, miss, dup, n_busco = parse_busco_summary(summary) if summary else (0.0, 0.0, 0.0, 0.0, 0)

        # 3. Clean Labels
        label = step
        label = label.replace("01_Unicycler_Raw", "Unicycler")
        label = label.replace("03_Medaka", "Medaka")
        label = label.replace("04_Polypolish", "Polypolish")
        label = label.replace("02_Racon_racon0", "Racon (Seed)")
        label = label.replace("02_Racon_racon1", "Racon (Rd 1)")
        label = label.replace("02_Racon_racon2", "Racon (Rd 2)")
        if "_" in label and "Racon" not in label: 
             label = label.split('_', 1)[-1].replace('_', ' ')

        rows.append({
            "Step": label,
            "N50_bp": n50,
            "BUSCO_C": comp
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT_TABLE, sep="\t", index=False)
    print(f"[+] Table saved: {OUT_TABLE}")

    if df.empty or df["BUSCO_C"].max() <= 0:
        print("[!] Skipping plot: No valid BUSCO data found.")
        sys.exit(0)

    # --- PROFESSIONAL PLOTTING ---
    sns.set_theme(style="whitegrid")
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Axis 1: N50 Structure (Line Plot - Blue)
    color_n50 = '#1f77b4'
    
    sns.lineplot(data=df, x='Step', y='N50_bp', marker='o', markersize=10, 
                 linewidth=3, color=color_n50, ax=ax1, label='N50 (Contiguity)', 
                 zorder=10, legend=False)
    
    ax1.set_ylabel('Assembly N50', color=color_n50, fontweight='bold', fontsize=12)
    ax1.tick_params(axis='y', labelcolor=color_n50)
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(False)

    # FIX: Restore Standard Scientific Notation (1e6)
    # This allows Matplotlib to handle micro-variations correctly
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    # Add Headroom to N50 Axis
    n50_min = df['N50_bp'].min()
    n50_max = df['N50_bp'].max()
    n50_span = n50_max - n50_min if n50_max != n50_min else n50_max * 0.1
    ax1.set_ylim(n50_min - (n50_span * 0.1), n50_max + (n50_span * 0.35))

    # Axis 2: BUSCO Biology (Bar Plot - Green)
    ax2 = ax1.twinx()
    color_busco = '#2ca02c'
    ax2.bar(df['Step'], df['BUSCO_C'], color=color_busco, alpha=0.3, width=0.5, label='BUSCO Completeness', zorder=1)
    
    ax2.set_ylabel('BUSCO Completeness (%)', color=color_busco, fontweight='bold', fontsize=12)
    ax2.tick_params(axis='y', labelcolor=color_busco)

    # Extend BUSCO scale for headroom
    if df['BUSCO_C'].min() > 80:
        ax2.set_ylim(80, 115) 
    else:
        ax2.set_ylim(0, 115)

    # --- UNIFIED LEGEND ---
    lines = [plt.Line2D([0], [0], color=color_n50, linewidth=3, marker='o')]
    bars = [plt.Rectangle((0,0),1,1, color=color_busco, alpha=0.3)]
    
    ax2.legend(lines + bars, ['N50 (Contiguity)', 'BUSCO Completeness'], 
               loc='upper left', frameon=True, fancybox=True, framealpha=0.9)

    plt.title(f"Optimization Trajectory: {ISO}", fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(OUT_FIG, dpi=300)
    print(f"[+] Dashboard saved: {OUT_FIG}")

if __name__ == "__main__":
    main()
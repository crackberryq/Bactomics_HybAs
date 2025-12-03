# Bactomics HybAs v8.4-Lite  
**Engineering-Grade Hybrid Genome Assembly for Construction Biotechnology**

⚠️ **Technical Prerequisites: Intermediate Level**

This workflow runs in a **Linux-based environment** (Linux, macOS, or Windows via WSL2).  
While the workflow is fully automated, users should have a working understanding of:

- Linux command-line navigation  
- Conda environments  
- Basic troubleshooting of environment/installation issues  

Some Conda installation issues (dependency conflicts, environment solver failures, shell initialization) may occur and are **unrelated to HybAs**. These require basic technical familiarity to resolve.

---

## 🧱 What is Bactomics?

**Bactomics** is a bioinformatics platform designed specifically for:

- Geotechnical engineers  
- Construction biotechnologists  
- Researchers applying Microbially Induced Calcium Carbonate Precipitation (MICP)

Where traditional pipelines focus on biological discovery, **Bactomics focuses on engineering verification**, ensuring:

- Reproducibility  
- Standardization  
- Traceability  
- Lineage-aware QC  
- Engineering-ready outputs  

---

## 🎯 Mission

Enable non-specialists to:

- Validate isolate identity (genome-wide)  
- Confirm urease/MICP-related pathways  
- Apply standard QC regardless of sequencing provider  
- Generate reproducible results across labs and industrial settings  

---

## 📦 HybAs v8.4-Lite

**Hybrid Assembly & Verification for Illumina + Nanopore reads**

This workflow performs:

- Contamination detection (Kraken2)  
- Target-taxid whitelist filtering (optional)  
- Hybrid assembly (Unicycler)  
- Triple polishing (Racon → Medaka → Polypolish)  
- Structural and functional annotation (Prokka)  
- Full QC aggregation (MultiQC)

Biomedical modules (AMR, plasmids) are **removed** to keep the workflow focused on *industrial + engineering relevance*.

---

## ✅ Validated Configuration (Exact Versions)

HybAs v8.4-Lite was developed and tested under the following versions:

| Component | Version |
|----------|---------|
| **OS** | Ubuntu 22.04.5 LTS (on WSL2) |
| **Snakemake** | **9.10.1** |
| **Conda** | **25.9.1** |
| **Python** | **3.11.13** |
| **Conda channels (order)** | bioconda → conda-forge → defaults |

For reproducibility, users should match these versions as closely as possible.

---

## 🛠️ Step 0 — Platform Setup (WSL2/macOS/Linux) + Conda

### 🪟 Windows 10/11 — Install WSL2

Open **PowerShell (Admin)** and run:

```powershell
wsl --install
```

This installs:

- WSL2  
- Ubuntu (default)  
- Required kernel components  

Reboot if prompted.

Launch **Ubuntu** from the Start Menu and update packages:

```bash
sudo apt update && sudo apt upgrade -y
```

### 🍎 macOS (Intel & Apple Silicon)

Open **Terminal** and continue to the Miniconda installation step below.

### 🐧 Linux (native)

Open a terminal and continue to the Miniconda installation step below.

---

## 📦 Install Miniconda (Python 3.11)

Download Miniconda from:  
https://docs.conda.io/en/latest/miniconda.html

Run the installer (example for Linux):

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

Reload your shell:

```bash
source ~/.bashrc    # Linux/WSL
# or
source ~/.zshrc     # macOS
```

(Optional) Pin Conda to the validated version:

```bash
conda install conda=25.9.1
```

---

## 🐍 Step 1 — Create the HybAs Environment

Create a dedicated environment for the workflow:

```bash
conda create -n hybas \
  -c bioconda -c conda-forge -c defaults \
  snakemake=9.10.1 python=3.11

conda activate hybas
```

❗ **Do NOT install additional tools into this environment.**  
HybAs automatically creates separate per-rule environments for Unicycler, Medaka, Prokka, BUSCO, Racon, etc.

---

## 📥 Step 2 — Clone the Repository

```bash
git clone https://github.com/crackberryq/Bactomics_HybAs.git bactomics
cd bactomics
```

---

## 📂 Step 3 — Required Input Folder Structure

Your project should look like:

```text
base_dir/
└── isolate_name/
    └── raw/
        ├── illumina/
        │   ├── sample_R1.fastq.gz
        │   └── sample_R2.fastq.gz
        └── nanopore/
            ├── sample.fastq.gz
            └── ...
```

✔ Illumina: must contain `_R1` and `_R2`  
✔ Nanopore: any `.fastq` or `.fq` accepted  

---

## ⚙️ Step 4 — Configure `config.yaml`

Example:

```yaml
base_dir: /home/user/bactomics
isolate: isolate3

target_taxid: 400634
threads: 12

keep_percent: 95
racon_rounds: 2
medaka_model: ''

busco_lineage: bacteria_odb10
run_kraken: true
```

---

## 🦠 BUSCO Lineage Setup

To list available datasets:

```bash
busco --list-datasets
```

Download a lineage manually:

```bash
mkdir -p db/busco
cd db/busco
busco --download bacteria_odb10
```

Create a BUSCO config file:

```text
db/busco/config.ini
```

Contents:

```ini
[busco]
datasets_dir = /absolute/path/to/db/busco
```

Export before running:

```bash
export BUSCO_CONFIG_FILE=db/busco/config.ini
```

---

## 🧪 Step 5 — Optional Dry Run (Highly Recommended)

```bash
snakemake -s Snakefile --use-conda --cores 4 -n
```

This checks:

- Folder structure  
- Config correctness  
- Tool environments  

**Without** executing any jobs.

---

## 🏃 Step 6 — Run the Workflow

### Full Pipeline

```bash
snakemake -s Snakefile --use-conda -p --cores 12
```

### Generate ONLY final assembly

```bash
snakemake --use-conda -p --cores 12 isolate3/work/assembly.final.fasta
```

### MultiQC only

```bash
snakemake --use-conda -p --cores 1 isolate3/reports/multiqc/multiqc_report.html
```

---

## 📊 Summary of Workflow Steps

1. **Raw read QC & merging**  
2. **Cleaning**  
3. **Assembly**  
4. **Polishing**  
5. **Annotation & QC**

---

## 📁 Output Overview

| File / Folder | Description |
|---------------|-------------|
| `work/assembly.final.fasta` | Final polished genome |
| `annotation/<isolate>.gff` | Prokka structural annotation |
| `annotation/<isolate>.faa` | Predicted proteins |
| `reports/multiqc/multiqc_report.html` | Interactive QC report |
| `reports/busco/` | BUSCO completeness |
| `reports/kraken2/` | Kraken2 reports |
| `logs/` | Execution logs |

---

## 🧾 Academic Citations

If you use HybAs v8.4-Lite, please cite:

**Software / Workflow**  
Goldstein et al. *Bactomics HybAs: A modular workflow for hybrid genome assembly and taxon-aware quality control.* (Submitted, 2025)

**Application**  
Goldstein et al. *Genome-Resolved Study of Indigenous Lysinibacillus Bioprotectants.* (Under review, 2025)

---

## ⚖️ License

Released under the **MIT License**.  
See `LICENSE` for full details.

---

## 🧬 Bactomics HybAs v8.4-Lite

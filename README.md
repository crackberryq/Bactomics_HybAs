# Bactomics HybAs v8.4-Lite  
**Engineering-Grade Hybrid Genome Assembly for Construction Biotechnology**

⚠️ **Technical Prerequisites: Intermediate Level**  
This workflow runs in a **Linux-based environment** (Linux, macOS, or Windows via WSL2).  
Although heavily automated, users should understand:

- Linux command-line basics  
- Conda environment management  
- How to interpret and troubleshoot installation / dependency errors  

Some Conda issues (solver conflicts, channel ordering, shell initialization) are **system-level** and not caused by HybAs itself.

---

### 💬 Need Help?

If you are not fully comfortable with Linux, Conda, or Snakemake, you may **copy and paste this entire README into ChatGPT**.  
The model can assist you with:

- Step-by-step installation  
- WSL2/macOS/Linux setup  
- Configuring `config.yaml`  
- Verifying your folder structure  
- Troubleshooting errors or dependency conflicts  
- Running the workflow and interpreting results  

For an even smoother experience, you can use the dedicated Bactomics assistant available on the GPT Store:

👉 **https://chatgpt.com/g/g-692fd2ee755881919c653d1db9929f92-bactomics-hybas-v8-4-lite-bioinformatics-tutor**

This ensures that even first-time users can run **HybAs v8.4-Lite** reliably and confidently.

---

⚡ Auto-Installer (Optional)

Bactomics includes an automated installer for users who prefer a quick, guided setup.
This script will:

- Install Miniconda (if missing)

- Configure Conda with strict channel priority

- Create the bactomics environment with the validated versions

- Install Snakemake, Python, and Mamba

- Clone or update the Bactomics repository

- Run a Snakemake dry-run to verify installation

- Run the installer with:
  
```bash
curl -LO https://raw.githubusercontent.com/crackberryq/Bactomics_HybAs/main/autoinstall.sh
bash autoinstall.sh
```

⚠️ Important Warnings

The script modifies your Conda configuration, including channel order.

- Only run it on a system you fully control (Linux, macOS, or WSL2).

- If you use complex Conda/HPC setups, prefer manual installation.

- The script does NOT download BUSCO or Kraken2 databases — you must configure or point to them manually.

- If config.yaml is incomplete, the dry-run will warn you, but this does not affect installation success.

✔ Recommended When

- You want the fastest fully working installation

- You are new to Conda or Snakemake

- You are using a clean Ubuntu / WSL2 / macOS environment

- You want guaranteed reproducibility

❌ Not Recommended When

- You maintain multiple or custom Conda environments

- You rely on HPC module systems

- You want to control every dependency manually

📘 Full Installer Documentation

The full installer documentation is here:

➡️ [autoinstall.md](https://github.com/crackberryq/Bactomics_HybAs/blob/main/autoinstall.md)


If you need help installing or running Bactomics, you can:

Copy/paste this README into ChatGPT for guided assistance

Or use the dedicated Bactomics tutor:

👉 https://chatgpt.com/g/g-692fd2ee755881919c653d1db9929f92-bactomics-hybas-v8-4-lite-bioinformatics-tutor

---

## 🧱 What is Bactomics?

**Bactomics** is a modular genomic analysis suite purpose-built for:

- Geotechnical engineers  
- Construction biotechnologists  
- Researchers applying **Microbially Induced Calcium Carbonate Precipitation (MICP)**  

Where many pipelines prioritize *biological discovery*,  
**Bactomics focuses on engineering verification**, ensuring:

- Genome-wide identity validation  
- Confirmation of key metabolic pathways (e.g., urease operon)  
- Traceability and reproducibility  
- Standardized QC across laboratories and sequencing providers  

---

## 🎯 Mission

To provide **standardized, engineering-oriented pipelines** that enable non-specialists to:

1. **Verify bacterial identity** at genome scale  
2. **Confirm MICP-relevant pathways** and safety-relevant traits  
3. **Standardize QC** across bio-cementation and infrastructure projects  
4. **Ensure reproducibility** across sequencing vendors and platforms  

---

## 📦 Bactomics HybAs (v8.4-Lite)

**HybAs v8.4-Lite** is a targeted **hybrid assembly + taxon-aware QC** workflow for Illumina + Nanopore reads, orchestrated via Snakemake.

Biomedical modules (AMR, plasmids) have been intentionally removed to maintain **industrial and engineering focus**.

---

## 🚀 Key Features

### 🔧 Targeted Assembly & Taxon Validation
- Kraken2 contamination detection  
- Optional KrakenTools whitelisting using `target_taxid`  
- BUSCO lineage-aware completeness assessment  

### 🧬 Hybrid Assembly Engine
- **Unicycler** hybrid assembler (ONT structure + Illumina accuracy)  
- Optimized for circular microbial genomes  

### 🔄 Triple Polishing Chain
- **Racon** (ONT-based polishing)  
- **Medaka** (ONT consensus correction)  
- **Polypolish** (Illumina-based error correction)  

### 📊 Engineering-Ready Outputs
- MultiQC master report  
- QUAST and BUSCO metrics  
- SeqKit read and assembly statistics  
- Prokka structural and functional annotation  

### ⚡ Lightweight Industrial Architecture
- No AMR screening  
- No plasmid reconstruction  
- Reproducible, lean, and suitable for engineering workflows  

---

## ✅ Validated Configuration (Exact Versions)

HybAs v8.4-Lite was developed and tested under the following versions:

| Component | Version |
|----------|---------|
| **OS** | Ubuntu 22.04.5 LTS (via WSL2) |
| **Conda** | 25.9.1 |
| **Snakemake** | 9.10.1 |
| **Python** | 3.11.13 |
| **Conda channels (priority)** | conda-forge → bioconda → defaults |

For maximum reproducibility, users are encouraged to match these versions as closely as possible.

---

# 🛠️ Step 0 — Install Conda and Prepare a Clean Environment

Bactomics should be installed in a **clean, dedicated Conda environment** to avoid conflicts with pre-existing tools.

### 0.1 Install Conda

#### 🪟 Windows 10/11 (via WSL2)

1. Install WSL2 (Ubuntu) from the Microsoft Store:  
   <https://aka.ms/wslstore>  
2. Launch Ubuntu and install Miniconda:

```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

#### 🍎 macOS (Intel / Apple Silicon)

```bash
cd ~
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
source ~/.zshrc
```

> Use the x86_64 installer if running on Intel macOS.

#### 🐧 Linux (Ubuntu / Debian / RHEL / others)

```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

---

### 0.2 Configure Conda (Strict Channel Priority)

This is recommended for stable and reproducible installs:

```bash
conda config --set channel_priority strict
conda config --add channels conda-forge
conda config --add channels bioconda
```

Verify:

```bash
conda config --show channel_priority
conda config --show channels
```

Expected:

```text
channel_priority: strict
channels:
  - conda-forge
  - bioconda
  - defaults
```

---

### 0.3 Create the HybAs Execution Environment

Create a dedicated environment for running the Snakemake workflow:

```bash
conda create -n hybas python=3.11 snakemake=9.10.1 mamba -y
conda activate hybas
```

Optional: create an additional **test** environment that mirrors development conditions:

```bash
conda create -n hybas_test python=3.11 snakemake=9.10.1 mamba -y
conda activate hybas_test
```

Check versions:

```bash
python --version
snakemake --version
```

Expected:

```text
Python 3.11.x
9.10.1
```

> ❗ **Do not install unrelated tools into `hybas`**.  
> HybAs will automatically create separate per-rule Conda environments for tools such as Unicycler, BUSCO, Prokka, Medaka, Racon, etc.

---

# 📥 Step 1 — Clone the Repository

```bash
cd ~
git clone https://github.com/crackberryq/Bactomics_HybAs.git bactomics
cd bactomics
```

---

# 📂 Step 2 — Required Input Folder Structure

Your project directory (`base_dir`) should be structured as:

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

- ✔ Illumina must include paired-end files: `_R1` and `_R2`  
- ✔ Nanopore reads may be `.fastq` or `.fastq.gz` with any filename pattern  

---

# ⚙️ Step 3 — Configure `config.yaml`

Example configuration file:

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

Key fields:

- `base_dir` — Root directory containing isolate subfolders  
- `isolate` — Name of the folder for the current sample  
- `target_taxid` — NCBI taxonomy ID for KrakenTools whitelist (optional)  
- `threads` — Number of CPU threads to use  
- `keep_percent` — Fraction of best ONT reads retained by Filtlong  
- `racon_rounds` — Number of Racon polishing rounds  
- `medaka_model` — Medaka model name (leave empty for default)  
- `busco_lineage` — BUSCO lineage dataset name  
- `run_kraken` — Enable/disable Kraken2 steps  

---

# 📦 Step 4 — BUSCO Lineage Setup

List available datasets:

```bash
busco --list-datasets
```

Download a specific lineage (example: `bacteria_odb10`):

```bash
mkdir -p db/busco
cd db/busco
busco --download bacteria_odb10
```

Create a BUSCO configuration file:

```text
db/busco/config.ini
```

With contents:

```ini
[busco]
datasets_dir = /absolute/path/to/db/busco
```

Export the BUSCO config path before running:

```bash
export BUSCO_CONFIG_FILE=db/busco/config.ini
```

---

# 📝 Parameter Summary

| Key           | Description                 | Default           |
|---------------|-----------------------------|-------------------|
| `base_dir`    | Project directory           | `/home/user/bactomics` |
| `isolate`     | Sample folder name          | `isolate3`        |
| `target_taxid`| Kraken whitelist taxon ID   | `None`            |
| `threads`     | CPU threads                 | `12`              |
| `keep_percent`| Filtlong read retention (%) | `95`              |
| `racon_rounds`| Racon polishing iterations  | `2`               |
| `medaka_model`| Medaka model string         | `''`              |
| `busco_lineage`| BUSCO dataset              | `bacteria_odb10`  |
| `run_kraken`  | Run Kraken2 contamination steps | `true`       |

---

# 💻 System & Resource Requirements

### 🌐 Internet

Required for:

- Conda environment creation  
- BUSCO dataset download  
- Kraken2 database download  

### 💾 Disk Space (Approximate)

- Workflow code + Conda envs: **≈ 3 GB**  
- Kraken2 Standard DB: **≈ 60 GB**  
- BUSCO datasets (per lineage): **≈ 50 MB**  
- Output per isolate: **≈ 500 MB**  

### 🧠 Memory (RAM)

- Minimum: **16 GB**  
- Recommended: **32–64 GB** for large or complex datasets  

---

# 🧪 Optional Dry Run

Before running full analysis, perform a dry run to check configuration and paths:

```bash
snakemake -s Snakefile --use-conda --cores 4 -n -p
```

This will:

- Validate folder structure  
- Parse the Snakefile and `config.yaml`  
- Check that Conda environments can be created  

No jobs are actually executed in dry-run mode.

---

# 🏃 Step 5 — Run the Workflow

### Full Pipeline

```bash
snakemake -s Snakefile --use-conda -p --cores 12
```

### Final Assembly Only

```bash
snakemake --use-conda -p --cores 12 isolate3/work/assembly.final.fasta
```

### MultiQC Report Only

```bash
snakemake --use-conda -p --cores 1 isolate3/reports/multiqc/multiqc_report.html
```

> Adapt `isolate3` to match the `isolate` specified in your `config.yaml`.

---

# 📊 Workflow Summary (Conceptual)

1. **Raw QC & Merging**  
   - FastQC on Illumina  
   - NanoPlot on ONT  

2. **Contamination Checking**  
   - Kraken2 classification (Illumina + ONT)  
   - Optional KrakenTools whitelist filtering based on `target_taxid`  

3. **Read Cleaning & Filtering**  
   - fastp trimming of Illumina reads  
   - Filtlong filtering of ONT reads  

4. **Hybrid Assembly**  
   - Unicycler hybrid assembly (ONT + Illumina)  

5. **Polishing**  
   - Racon (ONT)  
   - Medaka (ONT consensus)  
   - Polypolish (Illumina refinement)  

6. **QC & Annotation**  
   - QUAST assembly metrics  
   - BUSCO completeness  
   - SeqKit statistics  
   - Coverage estimation  
   - MultiQC aggregation  
   - Prokka genome annotation  

---

# 📁 Output Overview

Typical key outputs include:

| Path | Description |
|------|-------------|
| `work/assembly.final.fasta` | Final polished genome assembly |
| `annotation/<isolate>.gff`  | Genome annotation (GFF) |
| `annotation/<isolate>.faa`  | Predicted proteins (FASTA) |
| `reports/multiqc/multiqc_report.html` | Integrated QC report |
| `reports/busco/`           | BUSCO completeness results |
| `reports/kraken2/`         | Kraken2 classification reports |
| `reports/nanoplot_raw/` / `reports/nanoplot_filt/` | ONT QC plots and tables |
| `reports/seqkit/`          | Read and assembly statistics |
| `reports/coverage/`        | Approximate ONT coverage estimates |
| `logs/`                    | Snakemake and tool log files |

---

# 🧾 Academic Citations

If you use **Bactomics HybAs v8.4-Lite**, please cite:

**Software / Workflow**  
Goldstein et al. *Bactomics HybAs: A modular workflow for hybrid genome assembly and taxon-aware quality control.* (Submitted, 2025)

**Application**  
Goldstein et al. *Genome-Resolved Study of Indigenous Lysinibacillus Bioprotectants.* (Under review, 2025)

---

# ⚖️ License

Released under the **MIT License**.  
See `LICENSE` for full terms.

---

# 🧬 Bactomics HybAs v8.4-Lite

Hybrid assembly tailored for construction biotechnology and MICP engineering applications.

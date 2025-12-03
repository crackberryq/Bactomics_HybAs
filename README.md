# Bactomics HybAs v8.4‑Lite  
**Engineering‑Grade Hybrid Genome Assembly for Construction Biotechnology**

⚠️ **Technical Prerequisites: Intermediate Level**  
This workflow runs in a **Linux-based environment** (Linux, macOS, or Windows via WSL2).  
While fully automated, users should have working knowledge of:

- Linux command-line navigation  
- Conda environment management  
- Interpreting installation or dependency errors  

Some Conda-related system issues (solver conflicts, shell initialization errors) may occur and are **unrelated to HybAs**.  
These require intermediate technical troubleshooting.

---

# 🧱 What is Bactomics?

**Bactomics** is a modular bioinformatics suite purpose‑built for:

- Geotechnical engineers  
- Construction biotechnologists  
- Researchers working with **Microbially Induced Calcium Carbonate Precipitation (MICP)**  

Where traditional pipelines focus on *biological discovery*,  
**Bactomics focuses on engineering verification**—ensuring:

- Identity verification  
- Pathway confirmation  
- Traceability  
- Reproducibility across labs  
- Standardization for bio‑cementation workflows  

---

# 🎯 Mission

To provide **standardized, engineering‑oriented pipelines** enabling non‑specialists to:

1. **Verify** bacterial isolate identity with genome-wide resolution  
2. **Confirm** urease operon and other MICP-relevant pathways  
3. **Standardize** QC across sequencing vendors  
4. **Ensure reproducibility** across industrial and academic workflows  

---

# 📦 Bactomics HybAs (v8.4‑Lite)

A targeted **hybrid assembly + QC verification** pipeline for Illumina + Nanopore reads, controlled via Snakemake.

Biomedical modules (AMR, plasmids) were intentionally removed for **industrial relevance**.

---

# 🚀 Key Features (Merged & Updated)

### 🔧 Targeted Assembly & Lineage Validation
- Kraken2 contamination detection  
- Optional taxon‑whitelisting using target_taxid  
- BUSCO lineage‑specific completeness validation  

### 🧬 Hybrid Assembly Engine
- **Unicycler** integrates ONT structure with Illumina accuracy  
- Optimized for circular microbial genomes  

### 🔄 Triple Polishing Chain
- **Racon** (ONT polishing)  
- **Medaka** (ONT consensus correction)  
- **Polypolish** (Illumina final correction)  

### 📊 Engineering‑Ready Outputs
- MultiQC summary report  
- Prokka annotation  
- BUSCO, QUAST, coverage checks  
- Comprehensive QC summaries  

### ⚡ Lean Industrial Architecture
- No AMR modules  
- No plasmid reconstruction  
- Lightweight & reproducible  

---

# ✅ Validated Configuration (Exact Versions)

| Component | Version |
|----------|---------|
| **OS** | Ubuntu 22.04.5 LTS (via WSL2) |
| **Conda** | 25.9.1 |
| **Snakemake** | 9.10.1 |
| **Python** | 3.11.13 |
| **Conda channels** | bioconda → conda-forge → defaults |

For maximum reproducibility, match these versions closely.

---

## ⚠️ STEP 0 — Install Conda + Create a Clean Bactomics Environment (Strongly Recommended)

Bactomics should always be installed in a **clean, isolated Conda environment**.  
This prevents dependency conflicts with any previously installed bioinformatics tools.

---

### 0.1 Install Conda

#### Windows 10/11
Bactomics requires **WSL2 (Ubuntu 22.04)**.

1. Install WSL2 from Microsoft Store:  
   https://aka.ms/wslstore
2. Inside WSL2, install Miniconda:

```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
macOS (Intel or Apple Silicon)
bash
Copy code
cd ~
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
source ~/.zshrc
Linux (Ubuntu / Debian / CentOS / RHEL)
bash
Copy code
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
0.2 Configure Conda strict channel priority
This is required for reproducible and conflict-free installation.

bash
Copy code
conda config --set channel_priority strict
conda config --add channels conda-forge
conda config --add channels bioconda
Verify:

bash
Copy code
conda config --show channel_priority
conda config --show channels
Expected output:

vbnet
Copy code
channel_priority: strict
channels:
  - conda-forge
  - bioconda
  - defaults
0.3 Create the Bactomics execution environment
Recommended environment (validated configuration):

bash
Copy code
conda create -n bactomics python=3.11 snakemake=9.10.1 mamba -y
conda activate bactomics
To replicate internal testing environment exactly:

bash
Copy code
conda create -n test python=3.11 snakemake=9.10.1 mamba -y
conda activate test
Verify:

bash
Copy code
python --version
snakemake --version
mamba --version
Expected versions:

nginx
Copy code
Python 3.11.x
Snakemake 9.10.1
Mamba 1.x.x
0.4 Clone the repository
bash
Copy code
cd ~
git clone https://github.com/crackberryq/bactomics_hybas.git bactomics
cd bactomics
0.5 Verify installation (dry-run)
This checks:

correct folder structure

required envs/*.yml files

Snakefile syntax

Conda env creation

Run:

bash
Copy code
snakemake --use-conda -n -p --cores 4
Successful output looks like:

arduino
Copy code
Building DAG of jobs...
Job counts...
This was a dry-run.
If you see no FileNotFoundError messages, installation is correct.

# 🐍 Step 1 — Create the HybAs Environment (Pinned)

```bash
conda create -n hybas \
  -c bioconda -c conda-forge -c defaults \
  snakemake=9.10.1 python=3.11

conda activate hybas
```

❗ **Do NOT install additional tools into this environment.**  
HybAs automatically generates isolated tool environments (Unicycler, BUSCO, Medaka, Prokka, Racon, etc.)

---

# 📥 Step 2 — Clone the Repository

```bash
git clone https://github.com/crackberryq/Bactomics_HybAs.git bactomics
cd bactomics
```

---

# 📂 Step 3 — Required Input Folder Structure

```
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

✔ Illumina must include `_R1` + `_R2`  
✔ Nanopore may be `.fastq` or `.fastq.gz`  

---

# ⚙️ Step 4 — Configuration: `config.yaml`

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

# 📦 BUSCO Lineage Setup

List available datasets:

```bash
busco --list-datasets
```

Download manually:

```bash
mkdir -p db/busco
cd db/busco
busco --download bacteria_odb10
```

BUSCO config file:

```
db/busco/config.ini
```

Contents:

```ini
[busco]
datasets_dir = /absolute/path/to/db/busco
```

Export:

```bash
export BUSCO_CONFIG_FILE=db/busco/config.ini
```

---

# 📝 Parameter Summary

| Key | Description | Default |
|-----|-------------|---------|
| base_dir | Project directory | /home/user/bactomics |
| isolate | Sample folder name | isolate3 |
| target_taxid | Kraken whitelist taxon | None |
| threads | CPU threads | 12 |
| keep_percent | Read retention | 95 |
| racon_rounds | Polishing iterations | 2 |
| medaka_model | Medaka model | '' |
| busco_lineage | BUSCO dataset | bacteria_odb10 |
| run_kraken | Run Kraken2 filter | true |

---

# 💻 System & Resource Requirements

### **Internet**
Required for:
- Conda environment creation  
- BUSCO database download  
- Kraken2 database download  

### **Disk Space**
- Workflow + Conda envs: **≈ 3 GB**  
- Kraken2 Standard DB: **≈ 60 GB**  
- BUSCO datasets: **≈ 50 MB**  
- Output per isolate: **≈ 500 MB**  

### **Memory**
- Minimum: **16 GB**  
- Recommended: **32–64 GB**  

---

# 🧪 Step 5 — Optional Dry Run

```bash
snakemake -s Snakefile --use-conda --cores 4 -n
```

---

# 🏃 Step 6 — Run the Workflow

### Full Pipeline

```bash
snakemake -s Snakefile --use-conda -p --cores 12
```

### Final Assembly Only

```bash
snakemake --use-conda -p --cores 12 isolate3/work/assembly.final.fasta
```

### MultiQC Only

```bash
snakemake --use-conda -p --cores 1 isolate3/reports/multiqc/multiqc_report.html
```

---

# 📊 Workflow Summary

1. QC & merging  
2. Kraken2 contamination check  
3. Optional taxon whitelisting  
4. fastp trimming  
5. Filtlong filtering  
6. Unicycler hybrid assembly  
7. Racon → Medaka → Polypolish  
8. BUSCO, QUAST, coverage  
9. MultiQC + Prokka annotation  

---

# 📁 Output Overview

| Path | Description |
|------|-------------|
| `work/assembly.final.fasta` | Final genome assembly |
| `annotation/<isolate>.gff` | Gene annotation |
| `annotation/<isolate>.faa` | Protein sequences |
| `reports/multiqc/multiqc_report.html` | QC overview |
| `reports/busco/` | BUSCO results |
| `reports/kraken2/` | Kraken2 contamination reports |
| `logs/` | Snakemake & tool logs |

---

# 🧾 Academic Citations

**Software / Workflow**  
Goldstein et al. *Bactomics HybAs: A modular workflow for hybrid genome assembly and taxon‑aware quality control.* (Submitted, 2025)

**Application**  
Goldstein et al. *Genome‑Resolved Study of Indigenous Lysinibacillus Bioprotectants.* (Under review, 2025)

---

# ⚖️ License
Released under the **MIT License**.

---

# 🧬 Bactomics HybAs v8.4‑Lite
Hybrid assembly tailored for construction biotechnology and MICP engineering.


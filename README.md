# Bactomics: The Engineering Genomics Suite

> **⚠️ Technical Prerequisites: Intermediate Level**  
> This workflow runs in a Linux environment (Ubuntu or WSL). While the bioinformatic steps are automated, **working knowledge of the command line and Conda/Mamba package management is required.**  
>
> *Note: Users may occasionally encounter system-specific Conda installation bugs (e.g., dependency conflicts, shell initialization) that are unrelated to the Bactomics code. Troubleshooting these environmental issues requires intermediate technical expertise.*

**Bactomics** is a modular bioinformatics platform purpose-built for **geotechnitional engineers**, **construction biotechnologists**, and researchers working with **Microbially Induced Calcite Precipitation (MICP)**.

Traditional bioinformatics tools focus on biological discovery.  
**Bactomics focuses on engineering verification.**  
It ensures that raw sequencing data produced by external providers meets the **quality, reproducibility, and traceability** required for infrastructure-related biotechnology.

---

## 🎯 Mission Statement

To provide **standardized pipelines** that allow non-specialists to:

1. **Verify** bacterial isolate identity using genome-wide analysis  
2. **Confirm** critical metabolic pathways (e.g., urease operon)  
3. **Standardize** QC across bio-cementation projects  
4. **Ensure reproducibility** across labs and industrial workflows  

---

# 📦 Bactomics HybAs (v1.0 – HybAs v8.4-lite)

### Targeted Hybrid Assembly & Verification Workflow (Illumina + Nanopore)

**HybAs** is a Snakemake-controlled hybrid assembly workflow that implements **lineage-aware quality control**, ensuring that the genome assembled corresponds to the target taxon identified via 16S rRNA or environmental expectations.

---

## 🚀 Key Features

### 🔧 Targeted Assembly
- Validates organism identity via lineage expectations  
- Kraken2 contamination filtering  
- BUSCO lineage-specific completeness  

### 🧬 Hybrid Assembly Engine
- Illumina (accuracy) + Nanopore (structure)  
- Assembled using **Unicycler** for circular genomes  

### 🔄 Triple Polishing
- **Racon** (Nanopore polishing)  
- **Medaka** (consensus correction)  
- **Polypolish** (Illumina error correction)  

### 📊 Engineering-Ready Outputs
- MultiQC master report  
- Prokka annotation  
- QC summaries  

### ⚡ Streamlined (“Lite”) Architecture
Removes biomedical modules (AMR, plasmids) for industrial relevance.

---

## ✅ Tested Configuration

- **OS:** Ubuntu 22.04.5 LTS (via WSL2)  
- **Manager:** Conda 25.9.1  
- **Workflow:** Snakemake 9.10.1  
- **Python:** 3.11.13  

---

# 🛠️ Installation

## 1. Clone the repository

```bash
git clone https://github.com/crackberryq/bactomics_hybas.git bactomics
cd bactomics
```

## 2. Install Snakemake & Mamba

```bash
conda install -c conda-forge -c bioconda snakemake mamba
```

---

# 💻 System & Resource Requirements

### **1. Internet Access**
Required for:
- Installing Conda environments  
- Downloading BUSCO datasets  
- Downloading Kraken2 DB  

### **2. Disk Space**
- Code + envs: **~3 GB**  
- Kraken2 Standard DB: **~60 GB**  
- BUSCO datasets: **~50 MB**  
- Output per isolate: **~500 MB**  

### **3. Memory (RAM)**
- Minimum: **16 GB**  
- Recommended: **32–64 GB**  

---

# 📂 Input Folder Structure

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

Illumina: `_R1` and `_R2` required  
Nanopore: any `.fastq` or `.fastq.gz`

---

# ⚙️ Configuration (`config.yaml`)

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

# 📦 BUSCO Lineage Database

Set `busco_lineage` based on 16S identification.

### Auto-download
```bash
busco --list-datasets
```

### Manual download
```bash
busco --download bacteria_odb10
```

### Local storage
```bash
mkdir -p db/busco/
cd db/busco/
wget https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.tar.gz
tar -xvf bacteria_odb10.tar.gz
```

Create:

```
db/busco/config.ini
```

With:

```
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
| isolate | Sample folder | isolate3 |
| target_taxid | KrakenTools whitelist | None |
| threads | CPU threads | 12 |
| keep_percent | ONT read retention | 95 |
| racon_rounds | Polishing rounds | 2 |
| medaka_model | Medaka model | '' |
| busco_lineage | BUSCO dataset | bacteria_odb10 |
| run_kraken | Enable Kraken2 | true |

---

# 🏃 Running the Pipeline

### Full pipeline
```bash
snakemake --use-conda -p --cores 12
```

### Final assembly only
```bash
snakemake --use-conda -p --cores 12 isolate3/work/assembly.final.fasta
```

### MultiQC report only
```bash
snakemake --use-conda -p --cores 1 isolate3/reports/multiqc/multiqc_report.html
```

---

# 📊 Workflow Summary

1. QC & Merging  
2. Kraken2 contamination check  
3. Optional whitelisting  
4. fastp trimming  
5. Filtlong filtering  
6. Unicycler assembly  
7. Racon, Medaka, Polypolish  
8. BUSCO, QUAST, Coverage  
9. MultiQC + Prokka  

---

# 📄 Output Files

| Path | Description |
|------|-------------|
| `work/assembly.final.fasta` | Final genome |
| `reports/multiqc/multiqc_report.html` | QC report |
| `annotation/<isolate>.gff` | Genome annotation |
| `reports/busco/` | BUSCO results |
| `logs/` | Execution logs |

---

# 📄 Academic Citations

**Software:**  
Goldstein et al. *Bactomics HybAs: A modular workflow for hybrid genome assembly and taxon-aware quality control.* (Submitted, 2025)

**Application:**  
Goldstein et al. *Genome-Resolved Study of Indigenous Lysinibacillus Bioprotectants.* (Under review, 2025)

---

# ⚖️ License  
Released under the **MIT License**.

---

**Bactomics HybAs v8.4-lite**  
Hybrid assembly tailored for construction biotechnology workflows.

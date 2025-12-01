> **⚠️ Technical Prerequisites: Intermediate Level**
> This workflow runs in a Linux environment (Ubuntu or WSL). While the bioinformatic steps are automated, **working knowledge of the command line and Conda/Mamba package management is required.**
>
> *Note: Users may occasionally encounter system-specific Conda installation bugs (e.g., dependency conflicts, shell initialization) that are unrelated to the Bactomics code. Troubleshooting these environmental issues requires intermediate technical expertise.*

# Bactomics: The Engineering Genomics Suite

**Bactomics** is a modular bioinformatics platform purpose-built for **geotechnical engineers**, **construction biotechnologists**, and researchers working with **Microbially Induced Calcite Precipitation (MICP)**.

Traditional bioinformatics tools focus on biological discovery.  
**Bactomics focuses on engineering verification.**  
It ensures that raw sequencing data produced by external providers meets the **quality, reproducibility, and traceability** required for infrastructure-related biotechnology.

---

## 🎯 Mission Statement

To provide **standardized, engineering-grade pipelines** that allow non-specialists to:

1. **Verify** bacterial isolate identity using genome-wide analysis.  
2. **Confirm** critical metabolic pathways (e.g., urease operon).  
3. **Standardize** QC across bio-cementation projects.  
4. **Ensure reproducibility** across labs and industrial workflows.

---

# 📦 Bactomics HybAs (v1.0 – HybAs v8.4-lite)

### Targeted Hybrid Assembly & Verification Workflow (Illumina + Nanopore)

**HybAs** is a Snakemake-controlled hybrid assembly workflow that implements **lineage-aware quality control**, ensuring that the genome assembled corresponds to the target taxon identified via 16S rRNA or environmental expectations.

---

## 🚀 Key Features

### 🔧 Targeted Assembly
- Uses 16S-based lineage to validate organism identity.
- Enforces taxon-aware QC: Kraken2 contamination filtering, BUSCO lineage-specific completeness.

### 🧬 Hybrid Assembly Engine
- Illumina (accuracy) + Nanopore (structure)
- Assembled using **Unicycler** for gap-free, circular genomes.

### 🔄 Triple Polishing
- **Racon** (ONT polishing)
- **Medaka** (consensus correction)
- **Polypolish** (Illumina error correction)

### 📊 Engineering-Ready Outputs
- MultiQC master report
- Prokka annotation (GFF, FAA, GBK)
- Comprehensive QC summaries

### ⚡ Streamlined (Lite) Architecture
Removes biomedical tools (AMR, plasmids) to focus on **industrial and engineering relevance**.

---

# 🛠️ Installation

## 1. Clone Repository

```bash
git clone https://github.com/crackberryq/bactomics_hybas.git
cd bactomics_hybas
```

## 2. Install Snakemake & Mamba

```bash
conda install -c conda-forge -c bioconda snakemake mamba
```

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

### Illumina Detection Rules
- Filenames must contain `_R1` and `_R2`.

### Nanopore Detection Rules
- Any `.fastq` or `.fastq.gz` is accepted.

---

# ⚙️ Configuration (`config.yaml`)

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

# 📦 BUSCO Lineage Database (Important)

BUSCO requires lineage datasets such as `bacteria_odb10`, `bacillales_odb10`, etc.  
Set `busco_lineage` according to your **16S identification**.

### Option A — Auto-download (easy)
BUSCO will automatically fetch lineages:
```bash
busco --list-datasets
```

### Option B — Manual Download
```bash
busco --download bacteria_odb10
```

### Option C — Store Lineages Inside Project
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

Add:
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
|------|------------|---------|
| base_dir | Root project directory | /home/user/project |
| isolate | Sample folder name | isolate3 |
| target_taxid | Whitelist TaxID for KrakenTools | None |
| threads | CPU threads | 12 |
| keep_percent | % ONT reads retained (Filtlong) | 95 |
| racon_rounds | Number of Racon polishing rounds | 2 |
| medaka_model | Medaka model ('' = auto) | '' |
| busco_lineage | BUSCO dataset | bacteria_odb10 |
| run_kraken | Enable Kraken2 | true |

---

# 🏃 Running the Pipeline

### Full Pipeline

```bash
snakemake --use-conda --cores 12
```

### Build only final assembly

```bash
snakemake --use-conda --cores 12 isolate3/work/assembly.final.fasta
```

### Regenerate MultiQC report

```bash
snakemake --use-conda --cores 1 isolate3/reports/multiqc/multiqc_report.html
```

---

# 📊 Workflow Summary

1. **QC & Merging** (FastQC, SeqKit, NanoPlot)  
2. **Kraken2 Contamination Profiling**  
3. **Whitelisting via KrakenTools** (optional)  
4. **fastp** trimming (Illumina)  
5. **Filtlong** ONT filtering  
6. **Hybrid Assembly** (Unicycler)  
7. **Polishing:**  
   - Racon (ONT)  
   - Medaka  
   - Polypolish (Illumina)  
8. **Validation:**  
   - BUSCO completeness  
   - QUAST metrics  
   - Coverage check  
9. **Reporting:**  
   - MultiQC  
   - Prokka annotation  

---

# 📄 Output Files

| Path | Description |
|------|-------------|
| `work/assembly.final.fasta` | Final polished genome |
| `reports/multiqc/multiqc_report.html` | Master QC report |
| `annotation/<isolate>.gff` | Annotated genome |
| `reports/busco/` | BUSCO completeness |
| `logs/` | Execution logs |

---

# ❓ Troubleshooting

### BUSCO HTML missing
**Cause:** incorrect lineage (e.g., using "auto").  
**Fix:** match BUSCO dataset to 16S lineage.

### KrakenTools errors
Ensure:
- `db/kraken2_std_db/` is present  
- qc environment contains `extract_kraken_reads.py`

### Polypolish failure
Occurs if BAM is corrupted; delete and rerun mapping.

---

# 📄 Academic Citations

### Software:
Goldstein et al.,  
**Bactomics HybAs: A targeted hybrid assembly workflow enabling genomic verification of ureolytic bacteria in geotechnical engineering.**  
*Submitted to MethodsX (2025).*

### Application:
Goldstein et al.,  
**Valorizing Food Waste for Scalable, Low-Carbon Concrete Protection: A Genome-Resolved Study of Novel Indigenous Lysinibacillus Bioprotectants.**  
*Under Review (2025).*

---

# ⚖️ License

Released under the **MIT License**.

---

**Bactomics HybAs v8.4-lite** – Engineering-grade hybrid assembly for construction biotechnology.

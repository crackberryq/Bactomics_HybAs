# Bactomics: The Engineering Genomics Suite

> **⚠️ Technical Prerequisites: Intermediate Level**
> This workflow runs in a Linux environment (Ubuntu or WSL). While the bioinformatic steps are automated, **working knowledge of the command line and Conda/Mamba package management is required.**
>
> *Note: Users may occasionally encounter system-specific Conda installation bugs (e.g., dependency conflicts, shell initialization) that are unrelated to the Bactomics code. Troubleshooting these environmental issues requires intermediate technical expertise.*

**Bactomics** is a modular bioinformatics platform purpose-built for **geotechnical engineers**, **construction biotechnologists**, and researchers working with **Microbially Induced Calcite Precipitation (MICP)**.

Traditional bioinformatics tools focus on biological discovery.  
**Bactomics focuses on engineering verification.** It ensures that raw sequencing data produced by external providers meets the **quality, reproducibility, and traceability** required for infrastructure-related biotechnology.

---

## 🎯 Mission Statement

To provide **standardized pipelines** that allow non-specialists to:

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

### ✅ Tested Configuration
This workflow was successfully validated on the following system:
* **OS:** Ubuntu 22.04.5 LTS (via WSL2)
* **Manager:** Conda 25.9.1
* **Workflow:** Snakemake v9.10.1
* **Language:** Python v3.11.13

---

# 🛠️ Installation

## 1. Clone Repository
Clone the repository into a folder named `bactomics` to match the default configuration.

```bash
git clone [https://github.com/crackberryq/bactomics_hybas.git](https://github.com/crackberryq/bactomics_hybas.git) bactomics
cd bactomics
2. Install Snakemake & MambaBashconda install -c conda-forge -c bioconda snakemake mamba
💻 System & Resource Requirements1. Internet AccessRequired. The pipeline must be able to connect to:Anaconda Cloud (to install software environments).NCBI/BUSCO servers (to download reference datasets if not pre-cached).2. Disk Space (Estimates)Code & Environments: ~3 GBKraken2 Standard Database: ~60 GB (Critical: Ensure you have space for this!)BUSCO Database: ~50 MB (Negligible)Project Output: ~500 MB per bacterial isolate (varies by sequencing depth).3. Memory (RAM)Minimum: 16 GB (may fail on large genomes).Recommended: 32 GB - 64 GB.Note: Hybrid assembly (Unicycler) and polishing (Medaka) are RAM-intensive.📂 Input Folder StructurePlaintextbase_dir/
└── isolate_name/
    └── raw/
        ├── illumina/
        │   ├── sample_R1.fastq.gz
        │   └── sample_R2.fastq.gz
        └── nanopore/
            ├── sample.fastq.gz
            └── ...
Illumina Detection RulesFilenames must contain _R1 and _R2.Nanopore Detection RulesAny .fastq or .fastq.gz is accepted.⚙️ Configuration (config.yaml)Example:YAMLbase_dir: /home/user/bactomics
isolate: isolate3

target_taxid: 400634
threads: 12

keep_percent: 95
racon_rounds: 2
medaka_model: ''

busco_lineage: bacteria_odb10
📦 BUSCO Lineage Database (Important)BUSCO requires lineage datasets such as bacteria_odb10, bacillales_odb10, etc.Set busco_lineage according to your 16S identification.Option A — Auto-download (easy)BUSCO will automatically fetch lineages:Bashbusco --list-datasets
Option B — Manual DownloadBashbusco --download bacteria_odb10
Option C — Store Lineages Inside ProjectBashmkdir -p db/busco/
cd db/busco/
wget [https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.tar.gz](https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.tar.gz)
tar -xvf bacteria_odb10.tar.gz
Create:db/busco/config.ini
Add:[busco]
datasets_dir = /absolute/path/to/db/busco
Export:Bashexport BUSCO_CONFIG_FILE=db/busco/config.ini
📝 Parameter SummaryKeyDescriptionDefaultbase_dirRoot project directory/home/user/bactomicsisolateSample folder nameisolate3kraken_dbPath to Kraken2 databasedb/kraken2_std_dbtarget_taxidWhitelist TaxID for KrakenToolsNonethreadsCPU threads12keep_percent% ONT reads retained (Filtlong)95racon_roundsNumber of Racon polishing rounds2medaka_modelMedaka model ('' = auto)''busco_lineageBUSCO datasetbacteria_odb10🏃 Running the PipelineFull PipelineBashsnakemake --use-conda -p --cores 12
Build only final assemblyBashsnakemake --use-conda -p --cores 12 isolate3/work/assembly.final.fasta
Regenerate MultiQC reportBashsnakemake --use-conda -p --cores 1 isolate3/reports/multiqc/multiqc_report.html
💡 Pro Tip: Seeing the DetailsUse the -p flag (--printshellcmds) to see the exact shell commands being executed.The pipeline uses tee to print tool logs to the screen in real-time. If you want it to run silently in the background, remove -p and redirect stderr (e.g., snakemake ... > run.log 2>&1).📊 Workflow SummaryQC & Merging (FastQC, SeqKit, NanoPlot)Kraken2 Contamination Profiling 3. Whitelisting via KrakenTools (optional)fastp trimming (Illumina)Filtlong ONT filteringHybrid Assembly (Unicycler)Polishing: - Racon (ONT)MedakaPolypolish (Illumina)Validation: - BUSCO completenessQUAST metricsCoverage checkReporting: - MultiQCProkka annotation📄 Output FilesPathDescriptionwork/assembly.final.fastaFinal polished genomereports/multiqc/multiqc_report.htmlMaster QC reportannotation/<isolate>.gffAnnotated genomereports/busco/BUSCO completenesslogs/Execution logs❓ TroubleshootingBUSCO HTML missingCause: incorrect lineage (e.g., using "auto").Fix: match BUSCO dataset to 16S lineage.KrakenTools errorsEnsure:db/kraken2_std_db/ is presentqc environment contains extract_kraken_reads.pyPolypolish failureOccurs if BAM is corrupted; delete and rerun mapping.📄 Academic CitationsSoftware:Goldstein et al.,Bactomics HybAs: A targeted hybrid assembly workflow enabling genomic verification of ureolytic bacteria in geotechnical engineering. Submitted to MethodsX (2025).Application:Goldstein et al.,Valorizing Food Waste for Scalable, Low-Carbon Concrete Protection: A Genome-Resolved Study of Novel Indigenous Lysinibacillus Bioprotectants. Under Review (2025).⚖️ LicenseReleased under the MIT License.Bactomics HybAs v8.4-lite – Hybrid assembly tailored for construction biotechnology workflows.

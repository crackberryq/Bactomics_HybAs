# Bactomics HybAs v9.0.0-Aware

**Engineering-grade hybrid bacterial genome assembly, taxon-aware QC, provenance capture, and downstream genome analysis**

HybAs Aware is an isolate-centered **Snakemake workflow package for bacterial hybrid assembly, taxon-aware QC, provenance capture, and downstream computational genome analysis**. In the current repository layout, isolate data live directly under `bactomics/`, while the workflow code, configuration, app, CI assets, and the downstream computational analysis workflow live under `bactomics/hybas/`.

---

## 🧱 What is Bactomics?

Bactomics is the broader project workspace in which **HybAs** operates as the hybrid assembly and analysis component. In the current repository layout:

* `bactomics/` is the top-level working area
* each isolate is a direct child of `bactomics/`
* `bactomics/hybas/` contains the HybAs workflow package
* `bactomics/hybas/validation/` contains a second Snakemake workflow for downstream computational genome analysis of completed HybAs results

The project is framed as an engineering-oriented, reproducibility-focused workflow package with emphasis on standardized QC, provenance, and genome-level inspection.

---

## 🎯 What HybAs does

The main workflow in `Hybas.smk` is designed for **batch processing across isolates** and supports:

* hybrid assembly
* ONT-only assembly
* Illumina-only assembly
* stage-resolved polishing
* Kraken-guided read filtering and decision recording
* provenance capture
* environment snapshot capture
* checksums and tool-version recording
* per-isolate QC, summary, and reporting outputs


---

## 🚀 Key capabilities

### Hybrid and single-platform assembly

The main HybAs workflow automatically resolves the run mode for each isolate based on discovered inputs:

* `hybrid`
* `ont_only`
* `illumina_only`

### Taxon-aware control logic

The provided config and Snakefile support Kraken-aware configuration, including decision modes such as:

* `off`
* `manual`
* `auto`
* `aware`

The current main config uses `kraken_decision_mode: aware`.

### Multi-stage polishing

The workflow supports these polishing modes:

* `none`
* `racon`
* `medaka`
* `tripolish`

`tripolish` is defined in the Snakefile as a staged cascade:

1. Racon for the configured number of rounds
2. Medaka
3. Polypolish

### Reporting and provenance

The repository targets structured per-isolate outputs for QC, summary reporting, checksums, workflow provenance, tool-version tracking, and MultiQC aggregation.

### Downstream genome analysis

The second Snakemake workflow performs post-assembly computational analysis across multiple phases, rather than simply checking whether an assembly completed.

---

## 👥 Intended audience

This repository is written primarily for:

* collaborators
* reviewers
* advanced external users

The documentation should read as a **software manual with publication-aware framing**. It is not intended to function as a paper surrogate, and it is not written as a beginner-first tutorial.

---

## ⚠️ Technical prerequisites

This repository is designed for users who are comfortable working in a shell-based environment.

The primary documented setup path is based on **Ubuntu 22.04.5 LTS**, plus a controller Conda environment that launches either the Snakemake workflows or the Streamlit app. The per-rule tool environments referenced by the Snakefiles are then created automatically by Snakemake.

Users should be comfortable with:

* command-line navigation
* Git clone and pull workflows
* Snakemake execution basics
* Conda environment management
* reading workflow logs and troubleshooting dependency problems
* editing YAML and TSV configuration files

Some environment or solver issues can arise from Conda, channels, or local shell setup rather than from HybAs itself.

---

## 💬 Support

The documentation is structured to support guided setup and troubleshooting, including:

* step-by-step installation
* creation of the base execution environment
* setup guidance for the supported OS baseline
* help editing `config.yaml`
* help checking folder structure and sample sheets
* help understanding runtime errors and dependency conflicts
* help interpreting key workflow outputs

For an even smoother experience, you can use the dedicated Bactomics assistant available on the GPT Store:

👉 https://chatgpt.com/g/g-692fd2ee755881919c653d1db9929f92-bactomics-hybas-v8-4-lite-bioinformatics-tutor

This ensures that even first-time users can run HybAs v9.0.0-Aware reliably and confidently.
---

## 🗂 Example Repository layout

```text
bactomics/
├── <isolate_1>/
│   ├── raw/
│   │   ├── illumina/
│   │   └── nanopore/
│   ├── work/
│   ├── reports/
│   ├── annotation/
│   ├── logs/
│   ├── metadata/
│   ├── benchmarks/
│   └── validation/
├── <isolate_2>/
│   └── ...
└── hybas/
    ├── Hybas.smk
    ├── config.yaml
    ├── samples.tsv
    ├── environment.yaml
    ├── .condarc
    ├── app/
    │   └── hybas_streamlit_app.py
    ├── ci/
    │   ├── config.ci.yaml
    │   └── samples.ci.tsv
    └── validation/
        ├── HybAs_validation.smk
        ├── validation_config.yaml
        ├── scripts/
        ├── envs/
        └── config/
```

### Layout model

* `bactomics/` is the top-level working area
* each isolate is expected to be a direct child of `bactomics/` for this version
* `bactomics/hybas/` contains the workflow package
* `environment.yaml` is intended for the initial controller environment
* `.condarc` should be used for Snakemake launches
* `bactomics/hybas/validation/` contains the Downstream Computational Genome Analysis Workflow

---

## 🧬 Workflow overview

### 1. Main HybAs assembly workflow

The main workflow is isolate-centric. For each isolate, it:

1. reads isolate metadata from `samples.tsv`
2. discovers input files under `raw/illumina/` and `raw/nanopore/`
3. determines the run mode automatically
4. merges raw inputs by platform
5. performs raw QC and read-level preprocessing
6. optionally performs Kraken-based profiling and read-cleaning decisions
7. assembles with Unicycler
8. performs stage-resolved polishing with one of the supported polishing policies
9. produces final assembly outputs
10. runs downstream QC and summary steps
11. writes provenance, checksums, tool versions, and run summaries

### 2. Downstream Computational Genome Analysis Workflow

The workflow under `validation/` is a separate Snakemake DAG that operates on completed HybAs isolate outputs. It is phase-based:

* **Phase 1**: stage progression, structural summaries, BUSCO summaries, dnadiff comparisons, polishing-change summaries, retention and coverage summaries, and stage manifest generation
* **Phase 2**: stage-wise gene-space and BUSCO degradation or hotspot analysis
* **Phase 3**: genome-wide edit-density, hotspot ranking, spatial clustering, sequence-context, feature-overlap, repeat or rRNA or read-support analyses, and final interpretive summaries
* **Plotting and reporting**: per-isolate plots, dashboard-style figures, run summaries, MultiQC packaging, and report bundling

This second workflow is recommended as part of the standard full HybAs usage path, while remaining technically separable from the main assembly workflow.

### Phase purpose in plain language

* **Phase 1** tracks how the assembly changes across stages and summarizes structural and QC signals at each step.
* **Phase 2** examines how gene-space and BUSCO-related signals change across the assembly progression.
* **Phase 3** identifies where polishing edits cluster in the genome and adds biological context around those hotspots.

---

## 📥 Input requirements

### Isolate directory layout

Each isolate is expected under `bactomics/<isolate>/`.

The Snakefile resolves inputs from:

* `raw/illumina/`
* `raw/nanopore/`

### Illumina input discovery

Illumina files are discovered in `raw/illumina/` using patterns equivalent to:

* `*_R1*.fastq.gz`
* `*_R1*.fq.gz`
* `*_R1*.fastq`
* `*_R1*.fq`
* `*_R2*.fastq.gz`
* `*_R2*.fq.gz`
* `*_R2*.fastq`
* `*_R2*.fq`

The workflow expects matched R1 and R2 input sets for Illumina mode detection.

### ONT input discovery

ONT files are discovered in `raw/nanopore/` using patterns equivalent to:

* `*.fastq.gz`
* `*.fq.gz`
* `*.fastq`
* `*.fq`

### Sample sheet

The main sample sheet is a TSV file referenced by `config.yaml`. The code requires an `isolate` column and also accepts:

* `target_taxid`
* `busco_lineage`
* `medaka_model`

Example header:

```tsv
isolate	target_taxid	busco_lineage	medaka_model
```

Blank values are allowed for optional fields.

### Isolate naming rules

For reliable workflow behavior, isolate IDs should follow these rules:

* isolate IDs must be unique
* use only letters, numbers, underscores, and hyphens
* avoid spaces
* avoid shell-special characters
* the isolate folder name must exactly match the `isolate` value in `samples.tsv`

### Example isolate tree

```text
bactomics/
└── isolate_A01/
    └── raw/
        ├── illumina/
        │   ├── isolate_A01_R1.fastq.gz
        │   └── isolate_A01_R2.fastq.gz
        └── nanopore/
            └── isolate_A01_ont.fastq.gz
```

---

## ⚙️ Configuration

### Main workflow config

The main workflow reads `config.yaml`. Based on the provided file and Snakefile, the major configuration groups are:

#### Core execution

* `base_dir`
* `sample_sheet`
* `threads`
* `workflow_version`

#### Read preprocessing and polishing

* `keep_percent`
* `racon_rounds`
* `polish_mode`
* `medaka_model`

#### Taxonomy and lineage defaults

* `target_taxid`
* `busco_lineage`

#### Module toggles

* `run_kraken`
* `run_busco`
* `run_prokka`
* `run_multiqc`
* `run_polypolish`

#### Workflow gates

* `setup_only`
* `admin_mode`
* `ci_mode`
* `allow_ont_only`
* `allow_illumina_only`
* `require_target_taxid_for_kraken`

#### Database settings

* `db_root`
* `kraken_db_name`
* `busco_download_lineages`

#### Kraken-aware decision parameters

* `kraken_decision_mode`
* `kraken_auto_passthrough_min_percent`
* `kraken_auto_use_both_platforms`
* `kraken_include_children`
* `kraken_auto_prefer_species_min_percent`
* `kraken_aware_delta_min_percent`
* `kraken_aware_moderate_delta_min_percent`
* `kraken_aware_high_target_min_percent`
* `kraken_aware_moderate_target_min_percent`

### Controller environment config

The top-level `environment.yaml` is intended for the **initial execution environment** used to launch Snakemake or the Streamlit app.

Recommended controller environment identity:

* environment name: `hybas`
* scope: launch and control only, not full rule-tool dependency coverage

Recommended pinned controller-layer packages:

* `python=3.11`
* `snakemake=8.24.*`
* `streamlit=1.39.*`
* `pyyaml=6.0.*`
* `pandas=2.2.*`
* `tabulate=0.9.*`
* `rich=13.9.*`

Add `pip` only if the app or controller layer actually requires pip-installed extras.

This controller environment is distinct from the rule-specific environments under `envs/` and `validation/envs/`, which are expected to be created automatically during workflow execution.

### Downstream computational analysis config

The second workflow reads `validation/validation_config.yaml`. Important keys include:

* `base_dir`
* `sample_sheet`
* `threads`
* `validation_workflow_version`
* `validation_env`
* `multiqc_env`
* `phase1_script`
* `phase2_script`
* `phase3_script`
* `plot_script`
* `phase1_lineage`
* `phase1_threads`
* `phase2_window_bp`
* `phase3_window_bp`
* `phase3_top_n`
* `phase3_min_feature_overlap_bp`
* `phase3_barrnap_flank_bp`
* `phase3_motif_kmer`
* `phase3_motif_min_count`
* `require_gff`
* `enable_multiqc`
* `enable_per_isolate_multiqc`
* `enable_report_bundle`
* `multiqc_title`
* `multiqc_comment`
* `multiqc_config`
* `report_index_name`

### Notes on config behavior

The Streamlit app exposes structured controls for many, but not all, config values. In particular, not all Kraken-aware settings are surfaced in the current UI.

The validation config also contains execution-style toggles such as `run_phase1`, `run_phase2`, `run_phase3`, and `run_plots`. These should be documented carefully unless the corresponding Snakefile logic is confirmed to wire them directly into DAG selection.

---

## 💻 Installation model

The repository should be documented as having a **two-layer installation model**.

### Layer 1: base system and controller environment

This layer covers:

* the pinned Ubuntu baseline
* installation of Git and Conda
* cloning or pulling the repository
* creation of the initial Conda environment from `environment.yaml`
* activation of that environment before launching Snakemake or Streamlit

This controller environment is the environment that starts the workflow.

### Layer 2: automatic per-rule environments

The environments referenced by the main HybAs workflow and the Downstream Computational Genome Analysis Workflow are expected to be created automatically by Snakemake during execution.

That means the public docs should clearly distinguish between:

* the **bootstrap environment** used to run HybAs itself
* the **automatic rule environments** used by the workflow tools

### `.condarc` requirement

Snakemake launches should be documented with the repository `.condarc` in scope. The Streamlit app already injects `CONDARC` into workflow launches, and the CLI documentation should mirror that expectation.

---

## 🛠 Installation and first setup

### 1. Prepare the supported base system

The primary documented baseline is:

* **Ubuntu 22.04.5 LTS**

Additional environments such as macOS or Windows via WSL2 may be practical, but Ubuntu 22.04.5 LTS is the pinned reference platform for installation guidance.

### 2. Clone or update the repository

First-time setup should document cloning the repository. Existing installations should document pulling updates.

```bash
git clone <TBA-public-repo-url> bactomics
cd bactomics/hybas
```

For updates:

```bash
git pull
```

The final README should replace `<TBA-public-repo-url>` once the public repository URL is finalized.

### 3. Create the controller environment

The initial environment should be created from `environment.yaml`. This is the environment intended to launch:

* `snakemake`
* `streamlit`
* the HybAs control layer

```bash
conda env create -f environment.yaml
conda activate hybas
```

### 4. Use the repository `.condarc`

Snakemake launches should use the repository `.condarc`.

```bash
export CONDARC=$(pwd)/.condarc
```

The Streamlit app already mirrors this behavior when launching workflows.

### 5. Let Snakemake create the rule environments automatically

After the controller environment is active, the workflow-specific environments referenced by the Snakefiles are expected to be created automatically when running with `--use-conda`.

---

## 🛠 Running the main workflow

All example paths below are sanitized to relative paths for repository documentation.

### Recommended working directory

Run the main assembly workflow from inside `bactomics/hybas/`.

```bash
cd bactomics/hybas
export CONDARC=$(pwd)/.condarc
```

### Standard run

```bash
snakemake -s Hybas.smk --use-conda --cores 8
```

### Dry run

```bash
snakemake -s Hybas.smk -n --use-conda --cores 8
```

### Rerun incomplete jobs

```bash
snakemake -s Hybas.smk --use-conda --rerun-incomplete --cores 8
```

### Unlock a working directory

```bash
snakemake -s Hybas.smk --unlock
```

### CI profile

The repository contains a reduced CI configuration under `ci/`.

```bash
snakemake -s Hybas.smk \
  --configfile ci/config.ci.yaml \
  --use-conda \
  --cores 4
```

---

## 🧪 Running the Downstream Computational Genome Analysis Workflow

### Recommended working directory

The safest documented pattern is to run the second workflow from inside `bactomics/hybas/validation/`.

```bash
cd bactomics/hybas/validation
export CONDARC=$(pwd)/../.condarc
```

### Standard run

```bash
snakemake -s HybAs_validation.smk \
  --configfile validation_config.yaml \
  --use-conda \
  --cores 8
```

### Dry run

```bash
snakemake -s HybAs_validation.smk \
  -n \
  --configfile validation_config.yaml \
  --use-conda \
  --cores 8
```

### Rerun incomplete jobs

```bash
snakemake -s HybAs_validation.smk \
  --configfile validation_config.yaml \
  --use-conda \
  --rerun-incomplete \
  --cores 8
```

### Unlock the validation workflow directory

```bash
snakemake -s HybAs_validation.smk --unlock
```

### Running from `hybas/` instead of `hybas/validation/`

If you prefer to launch from `bactomics/hybas/`, use explicit relative paths to the validation Snakefile and config.

```bash
cd bactomics/hybas
export CONDARC=$(pwd)/.condarc
snakemake -s validation/HybAs_validation.smk \
  --configfile validation/validation_config.yaml \
  --use-conda \
  --cores 8
```

---

## 🖥 Streamlit control app

The repository includes a Streamlit application in `app/hybas_streamlit_app.py`.

From the provided code, the app currently supports:

* standard, admin, and CI operation modes
* resource-profile based thread scaling
* sample-sheet editing
* structured main-config editing
* structured downstream-analysis-config editing
* launch, stop, and unlock controls for both workflows
* runtime log display for both workflows

The app resolves the project root dynamically and expects the following workflow layout under the HybAs root:

* `Hybas.smk` or `HybAs.smk`
* `config.yaml`
* `samples.tsv`
* `ci/config.ci.yaml`
* `ci/samples.ci.tsv`
* `validation/HybAs_validation.smk`
* `validation/validation_config.yaml`

The app already uses the repository `.condarc` when launching Snakemake commands. In documentation, the app should be presented as a **workflow control interface** layered on top of the CLI.

### Launch command

Document the Streamlit app explicitly as:

```bash
streamlit run app/hybas_streamlit_app.py
```

---

## 📈 MultiQC integration

The provided MultiQC configuration indicates that HybAs exports a custom summary table into MultiQC. The custom content includes fields for:

* sample and run mode
* workflow version
* polishing mode
* target taxid and BUSCO lineage
* assembly size, contig count, N50, and GC
* BUSCO summary fields
* ONT retention and estimated coverage
* QC status fields
* tool-version summary columns

This means the main workflow is designed to write `metadata/multiqc_custom_content.tsv` per isolate and feed it into MultiQC reporting.

---

## 🔍 Provenance and reproducibility

The main Snakefile explicitly targets structured provenance capture. Based on the provided rules and target names, this includes:

* workflow provenance snapshots
* system information
* conda environment exports and package manifests
* database hash capture
* input and output checksums
* per-isolate tool-version summaries
* run summaries in both TSV and JSON forms

The downstream computational analysis workflow also targets:

* validation workflow provenance
* system information
* validation script checksums
* validation tool-version manifests
* run summaries in TSV and JSON

---

## ⚠️ Current repository notes

### Validation folder naming

The on-disk folder is named `validation`, but that naming should not be treated as a problem. It is simply the current folder name for the Downstream Computational Genome Analysis Workflow.

### Controller environment versus rule environments

The top-level `environment.yaml` is intended for the initial execution layer, not as the full dependency definition for every workflow tool. The rule-specific environments are separate and are expected to be created automatically by Snakemake.

### `.condarc` is part of the execution model

The repository `.condarc` should be part of documented Snakemake execution, both for CLI examples and for understanding how the Streamlit app launches workflows.

### Version strings are not yet fully harmonized

The provided files currently contain multiple version strings, including `HybAs-8.8-ieee-batch`, `HybAs-8.7-ci`, and a Streamlit banner showing `Version: 8.7`. These should be harmonized before release.

### Public repository URL and contact details remain pending

The final public GitHub URL and the public maintainer contact address are still `TBA` and should be filled in when the repository is prepared for release.

---

## ⚡ Quickstart

### Main workflow quickstart

```bash
cd bactomics/hybas
export CONDARC=$(pwd)/.condarc
snakemake -s Hybas.smk --use-conda --cores 8
```

### Downstream computational analysis quickstart

```bash
cd bactomics/hybas/validation
export CONDARC=$(pwd)/../.condarc
snakemake -s HybAs_validation.smk --configfile validation_config.yaml --use-conda --cores 8
```

### Streamlit quickstart

```bash
cd bactomics/hybas
conda activate hybas
streamlit run app/hybas_streamlit_app.py
```

---

## 📌 Project identity and governance

* **Repository name:** `HybAs-Aware`
* **Public-facing project name:** HybAs Aware
* **Short form in documentation:** HybAs
* **GitHub About tagline:** Engineering-grade hybrid bacterial genome assembly, taxon-aware QC, provenance capture, and downstream genome analysis.
* **License recommendation:** MIT
* **Citation file:** include `CITATION.cff`, pointing first to the software repository and later updated to the paper when public
* **Maintainer:** Jimi Goldstein
* **Affiliation:** Human Link
* **Public contact:** `TBA`
* **Contribution model:** controlled development with limited outside contributions

---


# USAGE.md

## HybAs Aware usage guide

This document explains how to **run HybAs Aware after installation**, how to choose the appropriate execution path, how to work with the sample sheet and isolate layout, and how to launch the downstream computational analysis workflow after main HybAs outputs are available.

This guide assumes that:

* the repository has already been cloned
* the controller environment has already been created from `environment.yaml`
* the controller environment can launch `snakemake` and `streamlit`
* the repository `.condarc` is available

For bootstrap setup, see `INSTALL.md`.

---

## 1. Usage model

HybAs Aware has **three practical entry routes**:

* CLI execution of the main workflow
* CLI execution of the Downstream Computational Genome Analysis Workflow
* Streamlit control app for launching and monitoring both workflows

The main workflow and the downstream analysis workflow are separate Snakemake DAGs. The second workflow is recommended as part of the standard full HybAs usage path, but it remains technically separable.

---

## 2. Repository working directories

### Main workflow working directory

The recommended directory for running the main HybAs workflow is:

```bash
bactomics/hybas/
```

### Downstream workflow working directory

The recommended directory for running the Downstream Computational Genome Analysis Workflow is:

```bash
bactomics/hybas/validation/
```

### Streamlit working directory

The recommended directory for launching the Streamlit app is:

```bash
bactomics/hybas/
```

---

## 3. Before running anything

Activate the controller environment and export the repository `.condarc`.

### Main workflow context

```bash
cd bactomics/hybas
conda activate hybas
export CONDARC=$(pwd)/.condarc
```

### Downstream workflow context

```bash
cd bactomics/hybas/validation
conda activate hybas
export CONDARC=$(pwd)/../.condarc
```

The `.condarc` export is part of the expected execution model and should not be skipped.

---

## 4. Input layout expected by the main workflow

Each isolate is expected under:

```text
bactomics/<isolate>/
```

Within each isolate, the workflow resolves inputs from:

* `raw/illumina/`
* `raw/nanopore/`

### Example isolate layout

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

### Illumina discovery rules

Illumina reads are discovered from `raw/illumina/` using R1/R2-aware patterns equivalent to:

* `*_R1*.fastq.gz`
* `*_R1*.fq.gz`
* `*_R1*.fastq`
* `*_R1*.fq`
* `*_R2*.fastq.gz`
* `*_R2*.fq.gz`
* `*_R2*.fastq`
* `*_R2*.fq`

### ONT discovery rules

ONT reads are discovered from `raw/nanopore/` using patterns equivalent to:

* `*.fastq.gz`
* `*.fq.gz`
* `*.fastq`
* `*.fq`

---

## 5. Sample sheet usage

The main workflow reads isolate-level metadata from `samples.tsv`, as specified by `config.yaml`.

The required column is:

* `isolate`

Supported optional columns include:

* `target_taxid`
* `busco_lineage`
* `medaka_model`

### Minimal example

```tsv
isolate	target_taxid	busco_lineage	medaka_model
isolate_A01		bacteria_odb10	
isolate_B02	1236	bacteria_odb10	
```

### Isolate naming rules

For reliable workflow behavior:

* isolate IDs must be unique
* use only letters, numbers, underscores, and hyphens
* avoid spaces
* avoid shell-special characters
* the isolate folder name must exactly match the `isolate` value in `samples.tsv`

If the isolate name in the sample sheet and the isolate directory name do not match exactly, the workflow should be expected to fail or misresolve paths.

---

## 6. How run mode is determined

The main workflow resolves the run mode automatically per isolate based on discovered input files.

Possible run modes are:

* `hybrid`
* `ont_only`
* `illumina_only`

### `hybrid`

Chosen when both paired Illumina inputs and ONT inputs are present.

### `ont_only`

Chosen when ONT inputs are present and paired Illumina inputs are absent.

### `illumina_only`

Chosen when paired Illumina inputs are present and ONT inputs are absent.

Acceptance of single-platform modes is controlled by:

* `allow_ont_only`
* `allow_illumina_only`

If an isolate has no valid reads under `raw/illumina/` and `raw/nanopore/`, the workflow should fail during input validation.

---

## 7. Main workflow execution patterns

### Standard full run

From inside `bactomics/hybas/`:

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

### CI mode

The repository includes a reduced CI configuration under `ci/`.

```bash
snakemake -s Hybas.smk \
  --configfile ci/config.ci.yaml \
  --use-conda \
  --cores 4
```

---

## 8. What the main workflow does during a run

For each isolate, the main workflow:

1. reads isolate metadata
2. discovers Illumina and ONT inputs
3. determines run mode
4. merges or prepares platform-specific input sets
5. performs QC and read preprocessing
6. optionally performs Kraken-based profiling and decision logic
7. assembles with Unicycler
8. polishes according to the configured polishing mode
9. writes final assembly outputs
10. generates reports, metadata, and provenance records

### Supported polishing modes

The main workflow supports:

* `none`
* `racon`
* `medaka`
* `tripolish`

`tripolish` is the staged sequence:

1. Racon
2. Medaka
3. Polypolish

---

## 9. Representative outputs from the main workflow

For each isolate, representative outputs include:

### Final assembly outputs

* `work/final/assembly.final.fasta`

### Intermediate assembly and polishing outputs

* `work/assembly/unicycler/`
* `work/assembly/stage_unicycler.fasta`
* `work/polishing/racon/`
* `work/polishing/medaka/consensus.fasta`
* `work/polishing/polypolish/stage_polypolish.fasta`

### Reports

* `reports/fastqc_pre/`
* `reports/nanoplot_raw/`
* `reports/nanoplot_filt/`
* `reports/kraken2/`
* `reports/fastp/fastp.html`
* `reports/fastp/fastp.json`
* `reports/quast/report.txt`
* `reports/quast/report.tsv`
* `reports/busco/`
* `reports/multiqc/multiqc_report.html`
* `reports/seqkit/`
* `reports/coverage/ONT_coverage.txt`

### Annotation outputs

When annotation is enabled:

* `annotation/<isolate>.gff`
* `annotation/<isolate>.faa`

### Metadata and provenance outputs

* `metadata/run_mode.txt`
* `metadata/run_summary.tsv`
* `metadata/run_summary.json`
* `metadata/qc_assessment.tsv`
* `metadata/multiqc_custom_content.tsv`
* `metadata/system_info.tsv`
* `metadata/workflow_provenance.tsv`
* `metadata/tool_versions.tsv`
* `metadata/checksums.tsv`
* `metadata/skip_reasons.tsv`
* `metadata/kraken_decision.tsv`

---

## 10. Configuration knobs that affect usage

The main workflow is controlled by `config.yaml`.

The most usage-relevant keys are:

### Core execution

* `base_dir`
* `sample_sheet`
* `threads`
* `workflow_version`

### Read preprocessing and polishing

* `keep_percent`
* `racon_rounds`
* `polish_mode`
* `medaka_model`

### Taxonomy and lineage defaults

* `target_taxid`
* `busco_lineage`

### Module toggles

* `run_kraken`
* `run_busco`
* `run_prokka`
* `run_multiqc`
* `run_polypolish`

### Workflow gates

* `setup_only`
* `admin_mode`
* `ci_mode`
* `allow_ont_only`
* `allow_illumina_only`
* `require_target_taxid_for_kraken`

### Kraken-aware settings

* `kraken_decision_mode`
* `kraken_auto_passthrough_min_percent`
* `kraken_auto_use_both_platforms`
* `kraken_include_children`
* `kraken_auto_prefer_species_min_percent`
* `kraken_aware_delta_min_percent`
* `kraken_aware_moderate_delta_min_percent`
* `kraken_aware_high_target_min_percent`
* `kraken_aware_moderate_target_min_percent`

If you need field-by-field descriptions, that material belongs in `CONFIGURATION.md`.

---

## 11. Streamlit app usage

The repository includes a Streamlit control interface at:

```text
app/hybas_streamlit_app.py
```

Launch it from inside `bactomics/hybas/` with the controller environment active:

```bash
conda activate hybas
streamlit run app/hybas_streamlit_app.py
```

### What the app is for

The app is a workflow control interface layered on top of the CLI. It is useful when you want:

* a structured launcher for standard, admin, and CI runs
* sample-sheet editing
* structured config editing
* validation-config editing
* workflow start, stop, and unlock controls
* runtime log viewing

### What the app is not for

The app does not replace:

* the need for a correct repository layout
* the need for a working controller environment
* the need for valid sample-sheet and input paths
* the need to understand basic Snakemake execution behavior

---

## 12. When to use CLI versus Streamlit

### Use the CLI when:

* you want exact command-line control
* you are debugging workflow behavior
* you are working in a remote shell session
* you need reproducible command history
* you are integrating HybAs into another scripted workflow

### Use the Streamlit app when:

* you want a control interface for launching and monitoring runs
* you want to edit sample sheets in a structured way
* you want a lighter operational interface for repeated runs

The CLI should still be treated as the canonical execution layer.

---

## 13. Running the Downstream Computational Genome Analysis Workflow

The second workflow is documented publicly as the **Downstream Computational Genome Analysis Workflow**.

It should generally be run **after** the main HybAs workflow has completed for the isolates of interest.

### Recommended working context

```bash
cd bactomics/hybas/validation
conda activate hybas
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

### Unlock the downstream workflow directory

```bash
snakemake -s HybAs_validation.smk --unlock
```

### Running it from `hybas/` instead

```bash
cd bactomics/hybas
conda activate hybas
export CONDARC=$(pwd)/.condarc
snakemake -s validation/HybAs_validation.smk \
  --configfile validation/validation_config.yaml \
  --use-conda \
  --cores 8
```

---

## 14. What the downstream workflow does

The downstream workflow is phase-based.

### Phase 1

Tracks how the assembly changes across stages and summarizes structural and QC signals at each step.

Representative outputs include:

* `validation/phase1/tables/stage_manifest.tsv`
* `validation/phase1/tables/busco_by_stage.tsv`
* `validation/phase1/tables/dnadiff_stage_vs_final.tsv`
* `validation/phase1/tables/dnadiff_prev_vs_next.tsv`
* `validation/phase1/metadata/phase1_metadata.tsv`

### Phase 2

Examines how gene-space and BUSCO-related signals change across the assembly progression.

Representative outputs include:

* `validation/phase2/tables/stage_gene_space_summary.tsv`
* `validation/phase2/tables/busco_hotspots.tsv`
* `validation/phase2/tables/busco_degradation_events.tsv`
* `validation/phase2/metadata/phase2_metadata.tsv`

### Phase 3

Identifies where polishing edits cluster in the genome and adds biological context around those hotspots.

Representative outputs include:

* `validation/phase3/tables/edit_summary_phase3.tsv`
* `validation/phase3/tables/window_density_phase3.tsv`
* `validation/phase3/tables/hotspots_phase3.tsv`
* `validation/phase3/tables/spatial_clustering_phase3.tsv`
* `validation/phase3/tables/classification_phase3.tsv`
* `validation/phase3/tables/final_summary_phase3.tsv`
* `validation/phase3/tables/final_summary_phase3.json`

### Plotting and reporting

Representative outputs include:

* `validation/plots/phase1/phase1_busco_trajectory.png`
* `validation/plots/phase2/phase2_cds_count.png`
* `validation/plots/phase3/phase3_spatial_vmr.png`
* `validation/plots/dashboards/dashboard_final_genome_integrity.png`
* `validation/multiqc/multiqc_report.html`
* `validation/report/`
* `validation/run_summary.tsv`
* `validation/run_summary.json`

---

## 15. When not to run the downstream workflow yet

Do not expect the downstream workflow to succeed if:

* the main HybAs workflow has not completed for the target isolate
* the expected final assembly is missing
* required upstream metadata or stage outputs are missing
* the validation config points to the wrong `base_dir` or sample sheet
* a required GFF is needed by your chosen analysis path and is not available

---

## 16. Typical usage sequences

### Sequence A: main workflow only

Use this when you only need assembly, QC, annotation, and standard reporting.

```bash
cd bactomics/hybas
conda activate hybas
export CONDARC=$(pwd)/.condarc
snakemake -s Hybas.smk --use-conda --cores 8
```

### Sequence B: full recommended path

Use this when you want the complete HybAs execution path including downstream computational analysis.

```bash
cd bactomics/hybas
conda activate hybas
export CONDARC=$(pwd)/.condarc
snakemake -s Hybas.smk --use-conda --cores 8

cd validation
export CONDARC=$(pwd)/../.condarc
snakemake -s HybAs_validation.smk --configfile validation_config.yaml --use-conda --cores 8
```

### Sequence C: app-driven usage

Use this when you want to operate through the Streamlit control interface.

```bash
cd bactomics/hybas
conda activate hybas
streamlit run app/hybas_streamlit_app.py
```

---

## 17. Common usage mistakes to avoid

### Forgetting to activate the controller environment

The CLI and Streamlit app both depend on the `hybas` controller environment.

### Forgetting to export `CONDARC`

The repository `.condarc` is part of the expected execution model for Snakemake runs.

### Running from the wrong directory

Main workflow commands and downstream workflow commands are easiest to manage when launched from their intended working directories.

### Mismatched isolate names

The isolate folder name and the `isolate` value in `samples.tsv` must match exactly.

### Treating the downstream workflow as independent of HybAs outputs

It is a separate DAG, but it still depends on outputs from the main workflow.

### Assuming the Streamlit app edits every possible config key

The current app exposes many settings, but not all possible config fields, especially not all Kraken-aware parameters.

---

## 18. Quick reference commands

### Main workflow dry run

```bash
cd bactomics/hybas
conda activate hybas
export CONDARC=$(pwd)/.condarc
snakemake -s Hybas.smk -n --use-conda --cores 4
```

### Main workflow full run

```bash
cd bactomics/hybas
conda activate hybas
export CONDARC=$(pwd)/.condarc
snakemake -s Hybas.smk --use-conda --cores 8
```

### Downstream workflow dry run

```bash
cd bactomics/hybas/validation
conda activate hybas
export CONDARC=$(pwd)/../.condarc
snakemake -s HybAs_validation.smk -n --configfile validation_config.yaml --use-conda --cores 4
```

### Downstream workflow full run

```bash
cd bactomics/hybas/validation
conda activate hybas
export CONDARC=$(pwd)/../.condarc
snakemake -s HybAs_validation.smk --configfile validation_config.yaml --use-conda --cores 8
```

### Streamlit launch

```bash
cd bactomics/hybas
conda activate hybas
streamlit run app/hybas_streamlit_app.py
```

---

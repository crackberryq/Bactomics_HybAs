## HybAs Aware troubleshooting guide

This document covers common failure modes and diagnostic patterns for **HybAs Aware** and its **Downstream Computational Genome Analysis Workflow**.

It is organized around the actual repository execution model:

* a controller environment created from `environment.yaml`
* Snakemake execution using the repository `.condarc`
* a main workflow launched from `bactomics/hybas/`
* a second workflow launched from `bactomics/hybas/validation/`
* isolate-centric input and output paths under `bactomics/<isolate>/`

This guide is meant to help you identify **where a problem originates** before you start changing code or configuration blindly.

For installation steps, see `INSTALL.md`. For run patterns, see `USAGE.md`. For field-level settings, see `CONFIGURATION.md`.

---

## 1. Start with the failure category

Before debugging, classify the problem into one of these buckets.

### A. Bootstrap or environment problem

Examples:

* `conda` not found
* `snakemake` not found
* `streamlit` not found
* controller environment fails to create
* solver conflicts during environment creation

### B. Repository layout or path problem

Examples:

* `Hybas.smk` not found
* `config.yaml` not found
* `samples.tsv` not found
* wrong `base_dir`
* isolate folder not found
* downstream workflow cannot find upstream outputs

### C. Input discovery problem

Examples:

* reads are present but the workflow says no valid inputs were found
* Illumina reads are not detected as pairs
* ONT reads are not discovered
* isolate appears in sample sheet but no data are found

### D. Workflow configuration problem

Examples:

* module unexpectedly skipped
* wrong run mode selected
* unsupported polishing behavior
* Kraken-related logic behaves unexpectedly
* downstream workflow uses the wrong sample sheet or wrong base directory

### E. Rule-level execution problem

Examples:

* a specific tool fails inside Snakemake
* Conda environment creation fails for one rule
* a report step fails while upstream assembly succeeds
* downstream phases fail after Phase 1 or Phase 2 completed

### F. Operational control problem

Examples:

* lock directory error
* interrupted run
* rerun behavior confusion
* Streamlit app launches but workflow does not start
* process stops but lock persists

If you identify the bucket first, troubleshooting becomes much faster.

---

## 2. Quick triage checklist

Run this checklist before making edits.

### Controller environment

* `conda activate hybas` works
* `python --version` works
* `snakemake --version` works
* `streamlit --version` works

### Repository context

From the main workflow directory:

```bash
cd bactomics/hybas
pwd
ls
```

You should see at least:

* `Hybas.smk`
* `config.yaml`
* `samples.tsv`
* `environment.yaml`
* `.condarc`
* `validation/`

### Conda execution model

```bash
export CONDARC=$(pwd)/.condarc
```

### First safe workflow check

```bash
snakemake -s Hybas.smk -n --use-conda --cores 4
```

If this fails, do not start with the Streamlit app. Fix the CLI path first.

---

## 3. Bootstrap and environment problems

### Problem: `conda: command not found`

Cause:

* Conda is not installed, or shell initialization is incomplete.

What to check:

```bash
which conda
conda --version
```

What to do:

* install Miniconda if missing
* initialize Conda for the current shell
* restart the shell session
* retry activation

### Problem: `snakemake: command not found`

Cause:

* the controller environment is not active
* `environment.yaml` has not been created yet
* the environment was created incorrectly

What to check:

```bash
conda info --envs
conda activate hybas
snakemake --version
```

What to do:

* activate `hybas`
* if that fails, recreate the controller environment from `environment.yaml`

### Problem: `streamlit: command not found`

Cause:

* controller environment is not active
* Streamlit is not included in the controller environment yet

What to check:

```bash
conda activate hybas
streamlit --version
```

What to do:

* verify that `environment.yaml` includes Streamlit
* recreate or update the controller environment if needed

### Problem: Conda solve fails or is extremely slow

Cause:

* channel priority is not configured consistently
* user-global Conda state conflicts with the repo’s expectations
* dependency pins are too restrictive or stale

What to check:

```bash
conda config --show channel_priority
conda config --show channels
```

Expected pattern:

* strict channel priority
* `conda-forge`
* `bioconda`
* `defaults`

What to do:

* restore strict channel priority
* retry environment creation
* keep the controller environment lean
* avoid manually installing unrelated tools into `hybas`

---

## 4. `.condarc` and Conda-resolution problems

### Problem: workflow works differently between CLI and app

Cause:

* CLI execution was started without the repository `.condarc`, while the Streamlit app injects `CONDARC` automatically.

What to check:

From `bactomics/hybas/`:

```bash
echo $CONDARC
```

What to do:

Use:

```bash
export CONDARC=$(pwd)/.condarc
```

From `bactomics/hybas/validation/` use:

```bash
export CONDARC=$(pwd)/../.condarc
```

### Problem: rule environments fail to resolve even though controller env works

Cause:

* controller environment is fine, but per-rule Conda resolution is failing
* rule env YAMLs may have stricter or conflicting pins

What to do:

* confirm `--use-conda` is being used
* confirm `CONDARC` is exported
* inspect the first failing rule and its env file
* do not assume `environment.yaml` covers the rule environments

---

## 5. Repository layout and path problems

### Problem: `Hybas.smk` not found

Cause:

* running from the wrong directory
* file name mismatch

What to check:

```bash
pwd
ls
```

Expected main workflow directory:

```text
bactomics/hybas/
```

What to do:

* move into `bactomics/hybas/`
* run `snakemake -s Hybas.smk ...`

### Problem: `HybAs_validation.smk` not found

Cause:

* running downstream commands from the wrong directory
* using the wrong relative path

What to do:

Either run from:

```bash
cd bactomics/hybas/validation
snakemake -s HybAs_validation.smk ...
```

or from:

```bash
cd bactomics/hybas
snakemake -s validation/HybAs_validation.smk ...
```

### Problem: sample sheet not found

Cause:

* `sample_sheet` path in config is wrong
* relative path is being interpreted from a different working context than expected

What to check:

* `config.yaml`
* `validation/validation_config.yaml`
* current working directory

What to do:

* confirm the path is correct relative to the Snakefile execution context
* keep published config paths relative and clean

### Problem: outputs appear under the wrong directory

Cause:

* `base_dir` is wrong
* sample sheet points at the wrong project root

What to check:

* `base_dir` in `config.yaml`
* `base_dir` in `validation/validation_config.yaml`

What to do:

* correct `base_dir`
* rerun a dry run before restarting the full workflow

---

## 6. Sample sheet and isolate problems

### Problem: isolate listed in sample sheet but workflow cannot find data

Cause:

* isolate folder name does not match the `isolate` value exactly
* isolate exists in the wrong parent directory
* `base_dir` is wrong

What to check:

* `samples.tsv`
* actual directory name under `bactomics/`

Rules to remember:

* isolate IDs must be unique
* avoid spaces and shell-special characters
* isolate folder name must exactly match the sample-sheet `isolate` value

### Problem: one isolate runs and another is skipped or fails unexpectedly

Cause:

* isolate-specific optional fields may differ
* one isolate has incomplete inputs
* run mode resolves differently between isolates

What to do:

* compare isolate directories directly
* compare the relevant sample sheet row values
* confirm that both isolates have the expected `raw/illumina/` and/or `raw/nanopore/` files

---

## 7. Input discovery problems

### Problem: workflow says no valid reads were found

Cause:

* files do not match the discovery patterns
* files are in the wrong directory
* Illumina reads are missing either R1 or R2

Expected directories:

* `raw/illumina/`
* `raw/nanopore/`

Expected Illumina-style patterns include:

* `*_R1*.fastq.gz`
* `*_R2*.fastq.gz`
* and equivalent `.fq` / uncompressed variants

Expected ONT-style patterns include:

* `*.fastq.gz`
* `*.fq.gz`
* `*.fastq`
* `*.fq`

What to do:

* verify the files are under the correct isolate directory
* verify Illumina files include both R1 and R2 matches
* rename files only if necessary and consistently

### Problem: ONT-only or Illumina-only isolate is rejected

Cause:

* single-platform mode exists, but the config may forbid it

What to check:

* `allow_ont_only`
* `allow_illumina_only`

What to do:

* enable the appropriate allowance if that execution path is intended

### Problem: wrong run mode appears in outputs

Cause:

* the workflow resolved the mode from discovered inputs, not from user expectation

What to do:

* inspect the actual files present in `raw/illumina/` and `raw/nanopore/`
* confirm paired Illumina reads exist as a valid pair
* confirm ONT files are not missing or misnamed

---

## 8. Main workflow configuration problems

### Problem: Kraken steps do not run

Cause:

* `run_kraken: false`
* taxon-related prerequisites not satisfied under current settings

What to check:

* `run_kraken`
* `target_taxid`
* `require_target_taxid_for_kraken`
* `kraken_decision_mode`

### Problem: BUSCO or Prokka steps do not run

Cause:

* corresponding module toggles are off

What to check:

* `run_busco`
* `run_prokka`

### Problem: MultiQC report is missing

Cause:

* `run_multiqc: false`
* upstream summary files are missing
* rule failure occurred before MultiQC aggregation

What to check:

* `run_multiqc`
* `metadata/multiqc_custom_content.tsv`
* upstream logs

### Problem: polishing behavior is not what you expected

Cause:

* `polish_mode` differs from your assumption
* `racon_rounds` is smaller or larger than intended
* Medaka model is blank or different than expected

What to check:

* `polish_mode`
* `racon_rounds`
* `medaka_model`
* `run_polypolish`

Remember:

* `tripolish` is a staged path, not a synonym for a single tool

---

## 9. Kraken-aware logic problems

### Problem: taxon-aware decision seems too permissive or too strict

Cause:

* the decision-policy thresholds are config-driven
* UI controls do not expose every advanced Kraken-aware field

What to check in `config.yaml`:

* `kraken_decision_mode`
* `kraken_auto_passthrough_min_percent`
* `kraken_auto_use_both_platforms`
* `kraken_include_children`
* `kraken_auto_prefer_species_min_percent`
* `kraken_aware_delta_min_percent`
* `kraken_aware_moderate_delta_min_percent`
* `kraken_aware_high_target_min_percent`
* `kraken_aware_moderate_target_min_percent`

What to do:

* inspect the file-based config directly
* do not assume the Streamlit UI exposes all decision-policy fields
* document threshold changes carefully for reproducibility

---

## 10. CI-profile confusion

### Problem: workflow behavior in CI mode does not match standard runs

Cause:

* `ci/config.ci.yaml` is intentionally reduced

Visible CI differences include:

* lower thread count
* fewer Racon rounds
* multiple optional modules disabled
* `ci_mode: true`

What to do:

* do not use CI outputs to judge full workflow behavior
* use the main config for normal analytical runs

---

## 11. Downstream workflow problems

### Problem: downstream workflow fails immediately

Cause:

* main HybAs outputs are missing
* `base_dir` or `sample_sheet` in `validation/validation_config.yaml` is wrong
* wrong working directory or wrong `.condarc` export

What to check:

* final assembly exists for the isolate
* sample sheet points at the intended isolates
* `base_dir` points at the real isolate root
* command is run from `bactomics/hybas/validation/` or with correct relative paths

### Problem: downstream workflow cannot find final assembly or annotation

Cause:

* upstream files were not produced
* filenames or directories differ from the expected structure
* GFF is required for the chosen path

What to check:

* `work/final/assembly.final.fasta`
* `annotation/<isolate>.gff`
* `require_gff`

### Problem: Phase 1 succeeds but later phases fail

Cause:

* downstream workflow is phase-based and later stages depend on Phase 1 tables, manifests, and metadata

What to do:

* inspect the first failing downstream rule
* inspect outputs under `validation/phase1/`
* do not debug Phase 3 first if Phase 1 outputs are already incomplete or malformed

### Problem: `run_phase1`, `run_phase2`, `run_phase3`, or `run_plots` do not behave as expected

Cause:

* these are config/UI-facing toggles, but you should not assume they fully prune the DAG unless the workflow logic explicitly confirms that behavior

What to do:

* treat these fields cautiously in troubleshooting
* inspect actual Snakemake target selection and rule graph behavior

---

## 12. Lock, interruption, and rerun problems

### Problem: Snakemake says the directory is locked

Cause:

* prior run exited unexpectedly
* workflow was interrupted
* app stopped without releasing lock cleanly

What to do for the main workflow:

```bash
cd bactomics/hybas
export CONDARC=$(pwd)/.condarc
snakemake -s Hybas.smk --unlock
```

What to do for the downstream workflow:

```bash
cd bactomics/hybas/validation
export CONDARC=$(pwd)/../.condarc
snakemake -s HybAs_validation.smk --unlock
```

### Problem: rerunning does not pick up incomplete jobs as expected

What to do:

Use:

```bash
--rerun-incomplete
```

This is the safest restart path when a run stopped mid-execution.

### Problem: you stopped a process from Streamlit but the directory still behaves as locked

Cause:

* process termination and lock cleanup are not always the same event

What to do:

* use the explicit unlock command after stopping the run

---

## 13. Streamlit app troubleshooting

### Problem: app opens but workflow does not start

Cause:

* controller environment is incomplete
* root layout is wrong
* paths expected by the app are missing
* Snakemake cannot be launched from the environment

What to check:

* `streamlit --version`
* `snakemake --version`
* presence of `Hybas.smk`, `config.yaml`, `samples.tsv`, `validation/HybAs_validation.smk`, and `validation/validation_config.yaml`

### Problem: app edits config but behavior still seems wrong

Cause:

* not every file-level setting is exposed in the UI
* advanced Kraken-aware fields remain file-based

What to do:

* inspect the YAML files directly
* treat file config as the source of truth

### Problem: sample sheet editor saves but run still uses wrong isolates

Cause:

* wrong operation mode selected
* CI mode may be using `ci/samples.ci.tsv`
* main mode uses `samples.tsv`

What to do:

* verify whether the app is in standard mode or CI mode
* confirm which sample sheet path is active

---

## 14. Output and reporting problems

### Problem: final assembly is missing

Cause:

* assembly or polishing failed upstream
* workflow terminated before finalization

What to inspect:

* `work/assembly/`
* `work/polishing/`
* isolate log files
* Snakemake terminal output

### Problem: MultiQC report is missing or incomplete

Cause:

* MultiQC module disabled
* upstream metadata table missing
* run failed before aggregation

What to inspect:

* `run_multiqc`
* `reports/multiqc/`
* `metadata/multiqc_custom_content.tsv`

### Problem: downstream plots are missing

Cause:

* plotting phase failed
* upstream phase outputs missing
* plotting toggle or report packaging path not aligned with the actual execution state

What to inspect:

* `validation/plots/`
* Phase 1, 2, and 3 tables
* first failing plotting rule

---

## 15. Provenance and version problems

### Problem: metadata show inconsistent version labels

Cause:

* the repository currently contains multiple version strings across configs and the Streamlit banner

What to do:

* harmonize `workflow_version`
* harmonize `validation_workflow_version`
* align the Streamlit display version before public release

### Problem: you cannot reproduce a prior run cleanly

Cause:

* config changed between runs
* sample sheet changed
* environment pins changed
* thresholds changed in file config but not recorded externally

What to do:

* preserve config snapshots per run
* preserve sample-sheet snapshots
* document changes to Kraken-aware thresholds and workflow version labels
* use metadata outputs such as provenance and tool-version tables for comparison

---

## 16. Safe diagnostic sequence

When a run fails, use this order.

### Step 1: confirm controller environment

```bash
conda activate hybas
python --version
snakemake --version
streamlit --version
```

### Step 2: confirm working directory

For the main workflow:

```bash
cd bactomics/hybas
pwd
ls
```

### Step 3: export `.condarc`

```bash
export CONDARC=$(pwd)/.condarc
```

### Step 4: dry run

```bash
snakemake -s Hybas.smk -n --use-conda --cores 4
```

### Step 5: inspect first real error only

Do not debug the 20th downstream symptom first. Start from the first failing rule or first missing file.

### Step 6: only then use unlock or rerun-incomplete

Use `--unlock` when lock state is the problem. Use `--rerun-incomplete` when interrupted jobs are the problem.

---

## 17. Common “don’t do this” list

* do not treat `environment.yaml` as the full tool stack for every rule
* do not skip exporting `CONDARC`
* do not assume Streamlit exposes every config key
* do not debug Phase 3 before checking whether Phase 1 outputs are valid
* do not use the CI config as if it represents the standard full workflow
* do not rename isolate directories without updating the sample sheet and related config
* do not assume a missing report means the reporting rule is broken; often the failure is upstream

---

## 18. Minimal recovery commands

### Main workflow dry run

```bash
cd bactomics/hybas
conda activate hybas
export CONDARC=$(pwd)/.condarc
snakemake -s Hybas.smk -n --use-conda --cores 4
```

### Main workflow unlock

```bash
cd bactomics/hybas
conda activate hybas
export CONDARC=$(pwd)/.condarc
snakemake -s Hybas.smk --unlock
```

### Main workflow rerun incomplete

```bash
cd bactomics/hybas
conda activate hybas
export CONDARC=$(pwd)/.condarc
snakemake -s Hybas.smk --use-conda --rerun-incomplete --cores 8
```

### Downstream workflow dry run

```bash
cd bactomics/hybas/validation
conda activate hybas
export CONDARC=$(pwd)/../.condarc
snakemake -s HybAs_validation.smk -n --configfile validation_config.yaml --use-conda --cores 4
```

### Downstream workflow unlock

```bash
cd bactomics/hybas/validation
conda activate hybas
export CONDARC=$(pwd)/../.condarc
snakemake -s HybAs_validation.smk --unlock
```

### Streamlit launch

```bash
cd bactomics/hybas
conda activate hybas
streamlit run app/hybas_streamlit_app.py
```

---



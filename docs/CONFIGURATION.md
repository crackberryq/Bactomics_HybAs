# CONFIGURATION.md

## HybAs Aware configuration reference

This document describes the configuration surfaces currently visible in the HybAs Aware repository.

It covers:

* the controller environment definition in `environment.yaml`
* the main workflow configuration in `config.yaml`
* the reduced CI configuration in `ci/config.ci.yaml`
* the downstream computational analysis configuration in `validation/validation_config.yaml`
* the relationship between file-based configuration and the Streamlit control interface

This file is intended as a **field-by-field configuration reference**, not as a full usage tutorial. For execution examples, see `USAGE.md`. For bootstrap setup, see `INSTALL.md`.

---

## 1. Configuration layers at a glance

The repository uses several distinct configuration layers.

### Controller layer

Defined in:

* `environment.yaml`
* `.condarc`

Purpose:

* defines the bootstrap environment used to launch Snakemake and Streamlit
* defines the Conda resolution behavior expected during workflow execution

### Main workflow layer

Defined in:

* `config.yaml`
* `samples.tsv`

Purpose:

* controls the main HybAs assembly and analysis workflow
* sets base paths, module toggles, polishing policy, taxon-aware behavior, and database locations

### CI layer

Defined in:

* `ci/config.ci.yaml`
* `ci/samples.ci.tsv`

Purpose:

* provides a reduced execution profile for CI or lightweight test runs

### Downstream analysis layer

Defined in:

* `validation/validation_config.yaml`

Purpose:

* controls the Downstream Computational Genome Analysis Workflow that runs on completed HybAs outputs

---

## 2. Controller-layer configuration

### `environment.yaml`

The top-level `environment.yaml` is intended for the **controller environment**, not as the full dependency specification for every workflow tool.

Its role is to define the environment that launches:

* `snakemake`
* `streamlit`
* controller-side Python dependencies

Recommended controller environment identity:

* environment name: `hybas`

Recommended controller-layer package set:

* `python=3.11`
* `snakemake=8.24.*`
* `streamlit=1.39.*`
* `pyyaml=6.0.*`
* `pandas=2.2.*`
* `tabulate=0.9.*`
* `rich=13.9.*`

Add `pip` only if the controller layer or app truly requires pip-installed extras.

### `.condarc`

The repository `.condarc` is part of the expected execution model.

It should be placed in scope for Snakemake launches using:

```bash
export CONDARC=$(pwd)/.condarc
```

or, from inside the downstream workflow directory:

```bash
export CONDARC=$(pwd)/../.condarc
```

This is important because the Streamlit app already injects `CONDARC` into workflow launches, so CLI usage should mirror that same behavior.

---

## 3. Main workflow configuration files

The main workflow is controlled primarily by:

* `config.yaml`
* `samples.tsv`

The workflow also discovers isolate data under the `base_dir` path.

---

## 4. `config.yaml` reference

The current main config is structured around several groups of settings.

### 4.1 Core execution fields

#### `base_dir`

Type: path-like string

Meaning:

* top-level working directory that contains isolate folders
* in the current repository model, isolates are direct children of `bactomics/`

Effect:

* all isolate path resolution is anchored under this directory
* incorrect values will cause sample folders and outputs to resolve incorrectly

#### `sample_sheet`

Type: path-like string

Meaning:

* path to the TSV file listing isolates and optional isolate-level overrides

Effect:

* determines which isolates are included in the workflow
* should point to a TSV with an `isolate` column

#### `threads`

Type: integer

Meaning:

* default core count used by the main workflow when no other CLI override is applied

Effect:

* influences workflow parallelism and resource allocation

#### `workflow_version`

Type: string

Meaning:

* human-readable workflow version label written into metadata and reporting outputs

Effect:

* useful for provenance and summary tables
* version strings should be harmonized across repo files before publication

---

### 4.2 Read preprocessing and polishing fields

#### `keep_percent`

Type: integer

Meaning:

* retention percentage used for ONT read filtering

Effect:

* affects how aggressively lower-priority ONT reads are discarded before downstream steps

#### `racon_rounds`

Type: integer

Meaning:

* number of Racon polishing rounds to perform when Racon is part of the selected polishing path

Effect:

* directly changes the number of polishing iterations

#### `polish_mode`

Type: string

Supported values confirmed by the workflow:

* `none`
* `racon`
* `medaka`
* `tripolish`

Meaning:

* selects the polishing strategy applied after assembly

Effect:

* `none`: no polishing stage beyond the initial assembly path
* `racon`: Racon-based polishing only
* `medaka`: Medaka-based polishing path
* `tripolish`: staged sequence of Racon, then Medaka, then Polypolish

#### `medaka_model`

Type: string

Meaning:

* optional Medaka model override

Effect:

* used when Medaka is part of the selected polishing path
* may be left blank if the workflow or environment uses default model behavior

---

### 4.3 Taxonomy and lineage defaults

#### `target_taxid`

Type: string or integer-like string

Meaning:

* default target taxonomy identifier used for taxon-aware or target-specific logic

Effect:

* relevant to Kraken-aware decision handling
* may be overridden at isolate level through the sample sheet if that pattern is used

#### `busco_lineage`

Type: string

Meaning:

* default BUSCO lineage dataset name

Effect:

* controls the lineage used for BUSCO-based completeness evaluation unless overridden

---

### 4.4 Module toggles

These fields switch major modules on or off in the main workflow.

#### `run_kraken`

Type: boolean

Meaning:

* enables or disables Kraken-related steps

#### `run_busco`

Type: boolean

Meaning:

* enables or disables BUSCO evaluation

#### `run_prokka`

Type: boolean

Meaning:

* enables or disables annotation via Prokka-related targets

#### `run_multiqc`

Type: boolean

Meaning:

* enables or disables MultiQC aggregation targets

#### `run_polypolish`

Type: boolean

Meaning:

* controls whether Polypolish-related steps are allowed where relevant to the selected polishing path

---

### 4.5 Workflow gates and mode controls

#### `setup_only`

Type: boolean

Meaning:

* used for setup-oriented workflow behavior rather than normal analytical execution

#### `admin_mode`

Type: boolean

Meaning:

* enables admin-oriented behavior and maintenance targets

#### `ci_mode`

Type: boolean

Meaning:

* flags the workflow as running in CI-style mode

#### `allow_ont_only`

Type: boolean

Meaning:

* controls whether ONT-only isolates are permitted to run

#### `allow_illumina_only`

Type: boolean

Meaning:

* controls whether Illumina-only isolates are permitted to run

#### `require_target_taxid_for_kraken`

Type: boolean

Meaning:

* controls whether Kraken-related execution requires a target taxid to be provided

---

### 4.6 Database settings

#### `db_root`

Type: path-like string

Meaning:

* root directory for workflow-managed database resources

Effect:

* used to resolve paths for Kraken, BUSCO, and other database-backed tooling as configured by the workflow

#### `kraken_db_name`

Type: string

Meaning:

* selects the Kraken database name under the configured database root

Effect:

* used by Kraken-related steps and integrity checks

#### `busco_download_lineages`

Type: list of strings

Meaning:

* lineages that may be referenced during BUSCO setup or maintenance operations

Effect:

* especially relevant in admin/setup contexts

---

### 4.7 Kraken-aware decision fields

These fields configure how taxon-aware logic behaves.

#### `kraken_decision_mode`

Type: string

Supported values confirmed by the workflow:

* `off`
* `manual`
* `auto`
* `aware`

Meaning:

* selects the decision strategy used for Kraken-aware behavior

The provided main config currently uses:

```yaml
kraken_decision_mode: aware
```

#### `kraken_auto_passthrough_min_percent`

Type: numeric

Meaning:

* threshold used by auto-style decision logic for direct acceptance or passthrough behavior

#### `kraken_auto_use_both_platforms`

Type: boolean

Meaning:

* controls whether evidence from both sequencing platforms is combined in auto-style logic

#### `kraken_include_children`

Type: boolean

Meaning:

* controls whether child taxa are included in target-aware interpretation

#### `kraken_auto_prefer_species_min_percent`

Type: numeric

Meaning:

* threshold used when auto logic prefers species-level interpretation

#### `kraken_aware_delta_min_percent`

Type: numeric

Meaning:

* higher-confidence delta threshold for aware-mode comparisons

#### `kraken_aware_moderate_delta_min_percent`

Type: numeric

Meaning:

* moderate-confidence delta threshold for aware-mode comparisons

#### `kraken_aware_high_target_min_percent`

Type: numeric

Meaning:

* high-confidence target-percentage threshold for aware-mode logic

#### `kraken_aware_moderate_target_min_percent`

Type: numeric

Meaning:

* moderate target-percentage threshold for aware-mode logic

These parameters should be documented as decision-policy controls rather than as universal biological truths.

---

## 5. `samples.tsv` reference

The sample sheet is a TSV file used by the main workflow.

### Required column

#### `isolate`

Type: string

Meaning:

* isolate identifier used to match the isolate directory under `base_dir`

Rules:

* must be unique
* should use only letters, numbers, underscores, and hyphens
* should not contain spaces or shell-special characters
* must exactly match the isolate folder name

### Supported optional columns

#### `target_taxid`

Type: string or integer-like string

Meaning:

* optional isolate-level override for taxon-aware logic

#### `busco_lineage`

Type: string

Meaning:

* optional isolate-level override for BUSCO lineage selection

#### `medaka_model`

Type: string

Meaning:

* optional isolate-level override for Medaka model choice

### Example

```tsv
isolate	target_taxid	busco_lineage	medaka_model
isolate_A01		bacteria_odb10	
isolate_B02	1236	bacteria_odb10	
```

---

## 6. CI configuration

The repository includes a reduced CI profile under:

* `ci/config.ci.yaml`
* `ci/samples.ci.tsv`

This is intended for reduced-footprint execution rather than full production-style runs.

### Key CI characteristics visible in the provided config

* `base_dir: ci`
* `sample_sheet: ci/samples.ci.tsv`
* `workflow_version: HybAs-8.7-ci`
* `threads: 4`
* `racon_rounds: 1`
* `run_kraken: false`
* `run_busco: false`
* `run_prokka: false`
* `run_multiqc: false`
* `run_polypolish: true`
* `ci_mode: true`

### Interpretation

This profile is a lighter-weight execution mode designed to test workflow behavior with fewer resources and fewer optional modules enabled.

It should not automatically be documented as equivalent to the standard full HybAs analytical path.

---

## 7. Downstream workflow configuration

The Downstream Computational Genome Analysis Workflow is configured through:

* `validation/validation_config.yaml`

This configuration controls the second Snakemake DAG that runs on completed HybAs outputs.

---

## 8. `validation/validation_config.yaml` reference

### 8.1 Core fields

#### `base_dir`

Type: path-like string

Meaning:

* root directory containing isolate folders and completed HybAs outputs

#### `sample_sheet`

Type: path-like string

Meaning:

* TSV describing the isolate set to be included in the downstream workflow

#### `threads`

Type: integer

Meaning:

* default core count for the downstream workflow

#### `validation_workflow_version`

Type: string

Meaning:

* version label written into downstream metadata and reports

---

### 8.2 Environment and script fields

#### `validation_env`

Type: path-like string

Meaning:

* path to the environment definition used by the downstream workflow core

#### `multiqc_env`

Type: path-like string

Meaning:

* path to the environment definition used for downstream MultiQC-related tasks

#### `phase1_script`

Type: path-like string

Meaning:

* entry script for Phase 1 processing

#### `phase2_script`

Type: path-like string

Meaning:

* entry script for Phase 2 processing

#### `phase3_script`

Type: path-like string

Meaning:

* entry script for Phase 3 processing

#### `plot_script`

Type: path-like string

Meaning:

* plotting entry script used for downstream figures and dashboards

---

### 8.3 Phase parameter fields

#### `phase1_lineage`

Type: string

Meaning:

* BUSCO lineage setting used by Phase 1 logic where applicable

#### `phase1_threads`

Type: integer

Meaning:

* phase-local thread setting for Phase 1 operations

#### `phase2_window_bp`

Type: integer

Meaning:

* window size in base pairs for Phase 2 window-based analyses

#### `phase3_window_bp`

Type: integer

Meaning:

* window size in base pairs for Phase 3 window-based analyses

#### `phase3_top_n`

Type: integer

Meaning:

* number of top-ranked regions or hotspots retained in key Phase 3 outputs

#### `phase3_min_feature_overlap_bp`

Type: integer

Meaning:

* minimum base-pair overlap used in Phase 3 feature-overlap logic

#### `phase3_barrnap_flank_bp`

Type: integer

Meaning:

* flank size used for Barrnap-associated context calculations

#### `phase3_motif_kmer`

Type: integer

Meaning:

* motif k-mer size used in Phase 3 motif-related calculations

#### `phase3_motif_min_count`

Type: integer

Meaning:

* minimum motif count used in Phase 3 motif-related logic

---

### 8.4 Execution and requirement fields

### Important behavior note

The config file contains execution-style toggles such as `run_phase1`, `run_phase2`, `run_phase3`, and `run_plots`. These should be documented carefully unless the workflow logic is fully confirmed to prune the DAG directly based on those fields.

#### `run_phase1`

Type: boolean

Meaning:

* intended as a Phase 1 execution toggle in config or UI

#### `run_phase2`

Type: boolean

Meaning:

* intended as a Phase 2 execution toggle in config or UI

#### `run_phase3`

Type: boolean

Meaning:

* intended as a Phase 3 execution toggle in config or UI

#### `run_plots`

Type: boolean

Meaning:

* intended as a plotting execution toggle in config or UI

#### `require_gff`

Type: boolean

Meaning:

* controls whether a GFF is required for the relevant downstream analysis path



---

### 8.5 Reporting and packaging fields

#### `enable_multiqc`

Type: boolean

Meaning:

* enables or disables downstream MultiQC packaging

#### `enable_per_isolate_multiqc`

Type: boolean

Meaning:

* enables or disables per-isolate MultiQC-related packaging behavior

#### `enable_report_bundle`

Type: boolean

Meaning:

* enables or disables downstream report bundle generation

#### `multiqc_title`

Type: string

Meaning:

* report title used for downstream MultiQC output

#### `multiqc_comment`

Type: string

Meaning:

* descriptive comment inserted into downstream MultiQC output

#### `multiqc_config`

Type: path-like string

Meaning:

* path to downstream MultiQC configuration content

#### `report_index_name`

Type: string

Meaning:

* file name used for the downstream report index

---

## 9. Streamlit app versus file-based configuration

The repository includes a Streamlit app at:

```text
app/hybas_streamlit_app.py
```

The app exposes structured controls for many configuration fields, but it is **not** a complete editor for every file-level parameter.

### Main workflow settings exposed in the app

The app exposes controls for many main config fields, including:

* `base_dir`
* `sample_sheet`
* `threads`
* `keep_percent`
* `racon_rounds`
* `busco_lineage`
* `medaka_model`
* `target_taxid`
* `run_kraken`
* `run_busco`
* `run_prokka`
* `run_multiqc`
* `run_polypolish`
* `setup_only`
* `admin_mode`
* `ci_mode`
* `allow_ont_only`
* `allow_illumina_only`
* `require_target_taxid_for_kraken`
* `polish_mode`
* `workflow_version`
* `db_root`
* `kraken_db_name`

### Important limitation

Not all Kraken-aware decision parameters are surfaced in the current UI. Advanced decision-policy tuning should still be treated as file-based configuration.

### Downstream settings exposed in the app

The app exposes several downstream fields, including:

* `isolate`
* `phase1_lineage`
* `phase2_window_bp`
* `phase3_window_bp`
* `phase3_top_n`
* `run_phase1`
* `run_phase2`
* `run_phase3`
* `run_plots`
* `enable_multiqc`

This is useful operationally, but the config file remains the source of truth for repository documentation.

---

## 10. MultiQC-related configuration

The repository includes a main MultiQC configuration file, currently represented by:

* `multiqc_config.yaml`

The provided MultiQC content indicates that HybAs writes a custom summary table into MultiQC using `metadata/multiqc_custom_content.tsv`.

Visible summary columns include fields such as:

* sample
* run mode
* workflow version
* polish mode
* target taxid
* BUSCO lineage
* assembly size
* contig count
* N50
* GC percent
* BUSCO summary values
* ONT retention and estimated coverage
* QC status fields
* tool-version summary fields

This means both config and metadata output naming should be kept stable if downstream reporting compatibility matters.

---

## 11. Configuration hygiene recommendations

### Keep paths relative in public documentation

The repository materials included absolute local paths. For publication-ready documentation, path examples should remain relative and sanitized.

### Keep isolate identifiers stable

Do not rename isolate directories after sample sheets and outputs have already been generated unless you also update the dependent config and metadata context.

### Harmonize version strings

The repository currently contains multiple version strings across files. These should be aligned before release.

### Treat controller and rule environments separately

Do not document `environment.yaml` as if it replaces the rule-specific environment files.

### Be cautious with UI assumptions

The Streamlit app is helpful operationally, but not every config parameter is visible there.

---

## 12. Common configuration mistakes to avoid

### Setting `base_dir` incorrectly

This can cause both input discovery and output resolution to fail.

### Mismatching `sample_sheet` and actual isolate folders

The isolate names in the sample sheet must correspond exactly to the isolate folder names.

### Assuming blank optional fields are always harmless

Blank values may be acceptable for some fields, but behavior depends on which modules are enabled and which execution path is chosen.

### Forgetting CI is a reduced profile

The CI config should not be assumed to represent full-featured execution.

### Assuming downstream phase toggles fully prune the DAG

Those settings should be documented carefully unless the pruning behavior is explicitly confirmed in workflow logic.

---

## 13. Minimal examples

### Example `config.yaml`

```yaml
base_dir: bactomics
sample_sheet: samples.tsv
threads: 8
keep_percent: 95
racon_rounds: 2
busco_lineage: bacteria_odb10
medaka_model: ""
target_taxid: ""
run_kraken: true
run_busco: true
run_prokka: true
run_multiqc: true
run_polypolish: true
setup_only: false
admin_mode: false
ci_mode: false
allow_ont_only: true
allow_illumina_only: true
require_target_taxid_for_kraken: false
polish_mode: tripolish
workflow_version: HybAs-Aware-TBA
db_root: db
kraken_db_name: kraken2_std_db
kraken_decision_mode: aware
```

### Example `samples.tsv`

```tsv
isolate	target_taxid	busco_lineage	medaka_model
isolate_A01		bacteria_odb10	
isolate_B02	1236	bacteria_odb10	
```

### Example `validation/validation_config.yaml`

```yaml
base_dir: bactomics
sample_sheet: samples.tsv
threads: 8
validation_workflow_version: HybAs-validation-TBA
phase1_lineage: bacteria_odb10
phase2_window_bp: 10000
phase3_window_bp: 10000
phase3_top_n: 20
run_phase1: true
run_phase2: true
run_phase3: true
run_plots: true
enable_multiqc: true
enable_report_bundle: true
```

---



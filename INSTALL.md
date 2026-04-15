## HybAs Aware installation guide

This document describes the **bootstrap installation path** for HybAs Aware.

It covers:

* the pinned reference operating system
* creation of the controller environment used to launch HybAs
* repository clone and update workflow
* use of the repository `.condarc`
* how automatic rule environments fit into the installation model
* first-launch checks for both CLI and Streamlit entry points

This guide does **not** attempt to replace workflow-specific usage documentation. Its purpose is to get a system into a state where HybAs can be launched correctly and reproducibly.

---

## 1. Installation model

HybAs Aware uses a **two-layer environment model**.

### Layer 1: controller environment

The top-level `environment.yaml` is intended for the controller layer. This is the environment that launches:

* `snakemake`
* `streamlit`
* the HybAs control layer

Recommended controller environment identity:

* environment name: `hybas`

Recommended controller-layer packages:

* `python=3.11`
* `snakemake=8.24.*`
* `streamlit=1.39.*`
* `pyyaml=6.0.*`
* `pandas=2.2.*`
* `tabulate=0.9.*`
* `rich=13.9.*`

This environment is **not** intended to contain every workflow tool used by HybAs.

### Layer 2: automatic rule environments

The tool-specific environments referenced by the Snakefiles under `envs/` and `validation/envs/` are expected to be created automatically by Snakemake when running with `--use-conda`.

That means the public install docs should distinguish clearly between:

* the **bootstrap environment** used to start HybAs
* the **automatic per-rule environments** used during execution

---

## 2. Supported installation baseline

The pinned reference platform for installation guidance is:

* **Ubuntu 22.04.5 LTS**

Additional environments such as macOS or Windows via WSL2 may be workable, but the main documented baseline is Ubuntu 22.04.5 LTS.

---

## 3. Before you start

This repository assumes you are comfortable with:

* shell navigation
* Git clone and pull workflows
* Conda environment activation
* Snakemake command execution
* editing YAML and TSV files when needed

This repository is aimed at collaborators, reviewers, and advanced external users. It is not written as a beginner-first installation tutorial.

---

## 4. Install the system prerequisites

On a clean Ubuntu 22.04.5 LTS system, ensure you have the basic command-line tools required to fetch and manage the repository.

A typical starting point is:

```bash
sudo apt update
sudo apt install -y git wget curl
```

If you already manage these tools another way, you do not need to reinstall them.

---

## 5. Install Conda

HybAs Aware expects a Conda-based controller environment.

If Conda is not already installed, install Miniconda first.

Example pattern:

```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

After installation, verify that Conda is available:

```bash
conda --version
```

If Conda is installed but not initialized in your shell, initialize it and restart the shell session before continuing.

---

## 6. Configure Conda channels

HybAs Aware uses Conda-managed environments, and stable resolution depends on predictable channel configuration.

Recommended channel setup:

```bash
conda config --set channel_priority strict
conda config --add channels conda-forge
conda config --add channels bioconda
```

Check the result:

```bash
conda config --show channel_priority
conda config --show channels
```

Expected pattern:

```text
channel_priority: strict
channels:
  - conda-forge
  - bioconda
  - defaults
```

The repository also includes a `.condarc` file. During workflow execution, Snakemake should be launched with the repository `.condarc` in scope.

---

## 7. Clone the repository

For first-time setup:

```bash
git clone <TBA-public-repo-url> bactomics
cd bactomics/hybas
```

### Updating an existing checkout

If the repository is already cloned locally:

```bash
cd bactomics/hybas
git pull
```

---

## 8. Inspect the top-level HybAs files

Before creating the environment, verify that the expected top-level workflow files are present inside `bactomics/hybas/`.

You should expect a structure like:

```text
Hybas.smk
config.yaml
samples.tsv
environment.yaml
.condarc
app/
ci/
validation/
```

The validation workflow is expected to remain under `validation/`. That folder name is simply the current repository convention.

---

## 9. Create the controller environment

The initial environment should be created from `environment.yaml`.

From inside `bactomics/hybas/`:

```bash
conda env create -f environment.yaml
conda activate hybas
```

This controller environment is intended to provide the launch layer for:

* Snakemake
* Streamlit
* controller-side Python dependencies

It is not the full dependency layer for every tool used by the pipeline.

### Recommended controller package set

The intended controller-layer package baseline is:

```text
python=3.11
snakemake=8.24.*
streamlit=1.39.*
pyyaml=6.0.*
pandas=2.2.*
tabulate=0.9.*
rich=13.9.*
```

Add `pip` only if the controller layer or app actually requires pip-installed extras.

---

## 10. Activate the repository `.condarc`

Snakemake runs should use the repository `.condarc` rather than relying only on the user-global Conda configuration.

From inside `bactomics/hybas/`:

```bash
export CONDARC=$(pwd)/.condarc
```

This is part of the actual execution model. The Streamlit app already mirrors this behavior when it launches workflow commands.

### Validation workflow path note

If launching from inside `bactomics/hybas/validation/`, use the parent `.condarc` path:

```bash
export CONDARC=$(pwd)/../.condarc
```

---

## 11. Verify the controller environment

After activating the environment, confirm that the core controller tools resolve correctly.

```bash
python --version
snakemake --version
streamlit --version
```

You should see:

* Python 3.11.x
* Snakemake 8.24.x
* Streamlit 1.39.x

Patch-level variation is acceptable if it matches the intended pinned series.

---

## 12. Understand what happens next during workflow execution

Once the controller environment is active and `CONDARC` is set, HybAs should be launched with `--use-conda`.

At that point, Snakemake is expected to create the rule-specific environments referenced by:

* the main HybAs workflow
* the Downstream Computational Genome Analysis Workflow

This means you do **not** need to manually preinstall every workflow tool into the `hybas` controller environment.

---

## 13. First dry-run check

Before attempting a full execution, run a dry run from inside `bactomics/hybas/`.

```bash
cd bactomics/hybas
export CONDARC=$(pwd)/.condarc
snakemake -s Hybas.smk -n --use-conda --cores 4
```

This is the safest first validation step for:

* the controller environment
* the Snakefile path
* the config load path
* the sample sheet path
* workflow graph construction

A dry run does not execute jobs, but it will still expose many configuration and path problems early.

---

## 14. First launch of the Streamlit control app

The repository includes a Streamlit control interface at:

```text
app/hybas_streamlit_app.py
```

From inside `bactomics/hybas/` with the controller environment active:

```bash
conda activate hybas
streamlit run app/hybas_streamlit_app.py
```

The app is a convenience control interface layered on top of the CLI. It does not replace the need for a working controller environment.

The current app supports:

* standard mode
* admin mode
* CI mode
* sample-sheet editing
* structured configuration editing
* workflow start and stop controls
* validation workflow start and stop controls
* runtime log viewing

---

## 15. Optional first launch of the downstream analysis workflow

The second Snakemake workflow is documented publicly as the **Downstream Computational Genome Analysis Workflow**.

The safest documented launch pattern is:

```bash
cd bactomics/hybas/validation
export CONDARC=$(pwd)/../.condarc
snakemake -s HybAs_validation.smk \
  --configfile validation_config.yaml \
  --use-conda \
  -n \
  --cores 4
```

This should usually be done after the main HybAs workflow has already produced the required isolate outputs.

---

## 16. Expected repository working model after installation

After installation, the intended working model is:

1. activate the `hybas` controller environment
2. move into `bactomics/hybas/`
3. export `CONDARC=$(pwd)/.condarc`
4. run HybAs with `--use-conda`
5. allow Snakemake to create the rule environments automatically
6. optionally launch the Streamlit control app from the same controller environment
7. run the downstream computational analysis workflow after main HybAs outputs are available

---

## 17. Common installation mistakes to avoid

### Treating `environment.yaml` as the full workflow dependency file

It is not. It should only define the controller environment used to launch HybAs.

### Forgetting to export `CONDARC`

The repository `.condarc` should be part of the documented execution flow.

### Installing all workflow tools manually into `hybas`

That defeats the intended per-rule environment model and makes troubleshooting harder.

### Running the downstream analysis workflow before main outputs exist

The second workflow depends on outputs from the main HybAs workflow.


---

## 18. Minimal installation checklist

Use this as a quick validation list.

* Ubuntu baseline prepared
* Git installed
* Conda installed and working
* channels configured
* repository cloned
* `bactomics/hybas/environment.yaml` present
* `bactomics/hybas/.condarc` present
* controller environment created
* `conda activate hybas` works
* `snakemake --version` works
* `streamlit --version` works
* `CONDARC` exported from the repository path
* `snakemake -s Hybas.smk -n --use-conda --cores 4` works

---


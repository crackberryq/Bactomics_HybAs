# Installation

Follow the [installation commands in the main README](../README.md#install). The controller environment is `hybas`, with Python 3.11 and Snakemake 9.6.2. Activate it before running the single `Hybas.smk` workflow.

Snakemake deploys and activates `envs/*.yml` with `--use-conda`. Post-deployment scripts check the required executables. There are no nested `conda run` calls and no pip-installed HybAs runner. Streamlit dependencies are optional and documented separately in the README.

Use `setup_envs --config setup_only=true` to prepare enabled environments without samples or databases. Normal runs require a populated sample sheet and all enabled databases. Database preparation is an explicit admin target; it does not occur during normal analysis.

Keep active files on a local Linux filesystem and configure the RAM budget available to Linux/WSL in `profiles/default/config.yaml`. Existing v8.8 outputs and environments should be preserved for comparison when adopting v9.

# Repository structure

| Path | Role |
|---|---|
| `Hybas.smk` | Directly runnable workflow and AWARE controller |
| `workflow_support.py` | Pure validation, stage planning and report helpers |
| `config.yaml`, `samples.tsv`, `samples.example.tsv` | Local run configuration and sample template |
| `environment.yaml`, `.condarc` | Controller installation and Conda channels |
| `envs/` | Tool environments and post-deployment executable checks |
| `profiles/default/config.yaml` | Scheduling resources and restart/log defaults |
| `report/` | Snakemake report captions |
| `multiqc_config.yaml` | Custom MultiQC table configuration |
| `scripts/check_database.py` | BUSCO dataset structure check used by admin preparation |
| `tests/` | Frozen-policy regressions and synthetic workflow contract tests |
| `.github/workflows/` | Regression and environment-deployment CI |
| `ci/` | Small example input profile |
| `app/` | Optional existing Streamlit interface |
| `docs/`, `CHANGELOG.md`, `CITATION.cff` | Documentation and release metadata |

Runtime results are written beneath the configured `base_dir`, not installed as Python package data. See the README for the per-isolate output layout. The separate historical `validation/` bundle is not included in this checkout.

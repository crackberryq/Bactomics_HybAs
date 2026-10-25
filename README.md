# HybAs 9.0.0-AWARE

[![Workflow regression](https://github.com/crackberryq/Bactomics_HybAs/actions/workflows/tests.yml/badge.svg?branch=HybAs_v9.0.0-AWARE)](https://github.com/crackberryq/Bactomics_HybAs/actions/workflows/tests.yml) [![Environment checks](https://github.com/crackberryq/Bactomics_HybAs/actions/workflows/environments.yml/badge.svg?branch=HybAs_v9.0.0-AWARE)](https://github.com/crackberryq/Bactomics_HybAs/actions/workflows/environments.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[Releases](https://github.com/crackberryq/Bactomics_HybAs/releases) · [Changelog](CHANGELOG.md) · [Security policy](SECURITY.md)

HybAs is a native Snakemake workflow for hybrid, Illumina-only and ONT-only bacterial assembly, with taxon-aware read handling, optional polishing, QC and provenance. Run `Hybas.smk` directly. Snakemake creates and activates the rule environments; there is no HybAs Python package to install.

This release fixes execution and reporting defects in the v8.8 source (`b14ff7688b72f1ae4a39bd99a4305d16a0e53ec1`). See [CHANGELOG.md](CHANGELOG.md) for behavior changes and validation limits. Preserve the old checkout and outputs when comparing versions.

## Install

Use Linux or WSL2, Bash and Conda. Keep the checkout, databases and active outputs on a local Linux filesystem. Synchronize completed results to OneDrive separately; cloud placeholders and concurrent synchronization can make active files unavailable.

```bash
git clone --branch HybAs_v9.0.0-AWARE https://github.com/crackberryq/Bactomics_HybAs.git HybAs
cd HybAs
export CONDARC="$PWD/.condarc"
conda env create -f environment.yaml
conda activate hybas
cp config.example.yaml config.yaml
cp examples/samples.example.tsv samples.tsv
# Edit config.yaml and samples.tsv for your workstation before normal analysis.
snakemake -s Hybas.smk setup_envs --use-conda --cores 4 --config setup_only=true
```

`environment.yaml` provides Python 3.11, Snakemake 9.6.2 and YAML support. Rule environments are in `workflow/envs/`. Their post-deployment checks fail if required executables are missing or resolve outside the managed environment. Optional Medaka, Polypolish and annotation environments are selected only when needed by the run (or requested during setup).

An existing environment can be updated with `conda env update -n hybas -f environment.yaml`. Keep your workstation configuration and sample sheet when upgrading; compare them with the new examples instead of overwriting them.

## Configure a run

Create your local `config.yaml` and `samples.tsv` from `config.example.yaml` and `examples/samples.example.tsv` (as shown under Install), then edit them locally. Both runtime files are gitignored to prevent accidental commits of workstation paths or sample metadata. Paths are resolved from the working directory. `base_dir` is the directory containing isolate folders; database paths are beneath `base_dir/db` unless configured otherwise.

Each sample requires an `isolate` column. Optional columns are `target_taxid`, `busco_lineage` and `medaka_model`. Isolate names must be unique and contain only letters, digits, dots, underscores or hyphens, beginning with a letter or digit.

Inputs live at:

- `<base_dir>/<isolate>/raw/illumina/*_R1*.fastq.gz` and matching `_R2` files. `.fq` and uncompressed FASTQ are also accepted. Lane names must match around the mate marker.
- `<base_dir>/<isolate>/raw/nanopore/*.fastq.gz` (also `.fq` and uncompressed FASTQ).

The workflow checks filenames, paired lanes, readability and the first FASTQ record before scheduling jobs. This is not a full read-pair or FASTQ integrity audit; decompression errors still fail the merge step.

A blank taxid or YAML `null` means no manual override. The literal strings `None` and `null`, booleans, zero and nonnumeric IDs are rejected. A populated sample-level ID overrides the global ID.

## AWARE policy

`kraken_decision_mode` accepts `off`, `manual`, `auto` and `aware`; it determines whether classification runs. `off` skips classification even if the legacy `run_kraken` flag is true.

The default AWARE policy is preserved:

1. Automatic election chooses the strongest inferred per-platform candidate.
2. Each platform's cleaning decision uses its own support for the elected target.
3. Cross-platform confidence remains conservative, including when one platform is absent. Low confidence does not itself override a platform's pass-through decision.

The confidence parameters now take effect: `kraken_aware_high_target_min_percent` (95), `kraken_aware_moderate_target_min_percent` (80), `kraken_aware_delta_min_percent` (20), and `kraken_aware_moderate_delta_min_percent` (10). Their defaults reproduce the previously active implementation. Changing them can change confidence labels; it does not change platform cleaning thresholds.

## Polishing routes

`racon_rounds` controls the number of Racon passes. A value of zero omits Racon.

| Requested mode | Hybrid | ONT only | Illumina only |
|---|---|---|---|
| `none` | Unicycler | Unicycler | Unicycler |
| `racon` | Racon | Racon | Unicycler |
| `medaka` | Racon → Medaka | Racon → Medaka | Unicycler |
| `tripolish` | Racon → Medaka → Polypolish | Racon → Medaka | Polypolish |

`run_polypolish: false` omits Polypolish and promotes the preceding stage. Each mate file is aligned independently with BWA's all-alignment option for Polypolish; secondary alignment records are retained. Paired BAMs used for mapping statistics are separate outputs.

If neither mate has mapped alignments, Polypolish records an explicit skip and retains its input assembly. A tool error fails the job. `metadata/polish_execution.json` records requested, planned, executed and skipped stages, stage hashes and the true final source stage. It does not infer execution from a configuration label.

## Resources

Edit `profiles/default/config.yaml` to fit the RAM available to Linux/WSL. Its initial total scheduling budget is **16000 MB**, with **one concurrent Kraken job**. This is a scheduling budget, not an operating-system memory limit. Leave RAM for the operating system and database/cache overhead.

- `assembly_mem_mb`: assembly job reservation, default 8192.
- `spades_memory_gb`: SPAdes memory argument forwarded by Unicycler, default 8; must fit the assembly reservation.
- `polishing_mem_mb`: polishing/mapping reservation, default 4096.
- `kraken_mem_mb`: optional explicit reservation. If omitted, HybAs estimates it from the Kraken core database files plus overhead. This estimate is not a measured peak.
- `min_fasta_length`: Unicycler export threshold, default 1000.

A large Kraken database can exceed the profile's example memory budget; size the budget and the database for your machine before running. HybAs rejects an assembly, Kraken or polishing reservation above that total before processing. Persistent Kraken processes/database reuse are not implemented in this release.

## Databases

Normal runs use existing local databases and fail early when required data are absent. Prepare them explicitly, after setting the paths and BUSCO lineages:

```bash
snakemake -s Hybas.smk download_busco_db --use-conda --cores 4 --config admin_mode=true
snakemake -s Hybas.smk download_kraken_db --use-conda --cores 8 --config admin_mode=true
```

Kraken preparation builds a standard database and can require substantial storage and memory. BUSCO preparation checks dataset structure before marking completion, downloads into a temporary location and retains an incomplete prior directory during replacement. A `.complete` file alone is insufficient evidence of a valid database.

## Run and resume

From the checkout with the intended local configuration:

```bash
conda activate hybas
export CONDARC="$PWD/.condarc"
snakemake -s Hybas.smk --use-conda --cores 8
```

Use `-n` for a dry run, `--printshellcmds` for command details, and `--keep-going` to let independent isolates finish after a failure. The default profile enables `--rerun-incomplete`. Detailed tool output is retained under each isolate's `logs/`; Snakemake reports the active stage and shows logs on failure.

After an interruption, first ensure the previous process has exited. If a stale lock remains:

```bash
snakemake -s Hybas.smk --unlock
snakemake -s Hybas.smk --use-conda --cores 8 --rerun-incomplete
```

Unlocking only releases a workflow lock. It does not repair an arbitrary broken environment or reset data. Snakemake rebuilds environments whose creation was incomplete; for a damaged environment marked complete, identify that specific environment from the failure before removing/recreating it. Do not delete the entire `.snakemake` directory or classification outputs to restart polishing.

## Outputs

Within each isolate directory:

- `work/final/assembly.final.fasta`: final promoted assembly.
- `work/assembly/` and `work/polishing/`: intermediate stages.
- `annotation/`: optional Prokka output.
- `reports/`: QUAST, BUSCO, read QC, coverage and MultiQC.
- `metadata/`: AWARE decisions, stage execution, hashes, run summaries and QC status.
- `logs/` and `benchmarks/`: diagnostics and resource measurements.

Global provenance is stored under `<base_dir>/provenance/`, including effective configuration, source hashes and deployed environment manifests. MultiQC explicitly depends on the custom summary and its configuration, so a summary refresh also refreshes the report.

```bash
snakemake -s Hybas.smk --report hybas_report.html
```

## Tests and qualification

```bash
python -m pip install pytest
python -m pytest tests -q
```

Controller tests compare default decisions and confidence with the frozen active v8.8 functions over 5,000 seeded cases, then exercise configuration and validation helpers. Workflow tests use real Snakemake with synthetic FASTQ and tool-contract stubs to check routing, separate SAM inputs, failure propagation, restart reuse, metadata and report refresh. These tests do **not** establish biological accuracy or equivalence of new assemblies to the published benchmark.

GitHub Actions runs the tests on Python 3.11 and 3.12 and checks Conda environment deployment. A hosted sandbox whose `/proc` exposes a different PID namespace may run `HYBAS_TEST_DISABLE_BENCHMARKS=1 python -m pytest tests -q`; that setting strips benchmark directives only from a temporary test copy. CI and the shipped workflow retain benchmarking.

The historical phase 1–3 scientific validation bundle is not present in this repository. It is separate from these software regression tests. The optional UI disables validation launch/save controls until that bundle is installed. Representative real-tool runs and comparisons to frozen biological outputs remain release-qualification work.

## Optional UI

The command-line workflow is the primary interface. To use the existing Streamlit dashboard, install its optional dependencies into `hybas`:

```bash
conda install -n hybas -c conda-forge streamlit pandas psutil
streamlit run app/hybas_streamlit_app.py
```

The dashboard launches Snakemake directly. It is not required for setup, execution or recovery.

## Citation and license

See [CITATION.cff](CITATION.cff) for citation metadata and [LICENSE](LICENSE) for license terms. Publication benchmark results belong to the version and inputs used to generate them; this update does not retroactively revalidate them.

## Repository maintenance and releases

The default branch may remain on a historical version until the repository owner changes it in **Settings → General → Default branch**. For v9.0, use the explicitly versioned branch shown in the clone command. The current CI checks validate software orchestration with mocks; they do not substitute for real-tool assembly qualification.

Repository contributors should use feature branches and pull requests; see [CONTRIBUTING](.github/CONTRIBUTING.md) and [SECURITY](SECURITY.md). Release notes for this version are in [docs/releases/v9.0.0-AWARE.md](docs/releases/v9.0.0-AWARE.md). To publish a GitHub Release, create a tag at the validated commit and publish it from the Releases page after CI passes.

## Experimental HybAs Research Report (standalone; no analysis rerun)

An optional offline, BaSIS-inspired HTML report can be generated from an **existing** isolate folder. This reporting-only tool does not modify the Snakefile, rerun Kraken, assemble genomes, or replace MultiQC.

```bash
conda activate hybas
python workflow/scripts/hybas_research_report.py \
  --base-dir /home/jimi/bactomics/TB_RUN \
  --isolate isolate001 \
  --audit-reads
```

Open `<base_dir>/<isolate>/reports/hybas_research/hybas_report.html`. The report shows AWARE target election, independent platform actions, confidence, assembly and polishing evidence, and an optional disk-backed read-ID audit. It reuses saved Kraken assignments and reconciles original merged FASTQs with cleaned FASTQs; it **does not classify reads again**.

Outputs include `aware_filter_audit.json`, `aware_removed_taxa.tsv`, and `aware_removed_categories.tsv` for Illumina, plus `aware_ont_filter_audit.json`, `aware_ont_removed_taxa.tsv`, and `aware_ont_removed_categories.tsv` for ONT when inputs are available. Human, bacterial, viral, fungal and other labels are **Kraken taxonomic assignments**, not confirmed biological origins. The audit counts read pairs, not independent mates. An integrity PASS verifies observed retention and assignment reconciliation, **not** biological correctness or every descendant whitelist decision.

For large datasets the temporary SQLite audit needs disk space and time. Omit `--audit-reads` for a quick evidence-only HTML report. This experimental tool is independent of the production DAG; do not audit files while an active job is writing them.

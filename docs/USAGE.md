# Usage

Edit `config.yaml` and `samples.tsv` in the checkout, then:

```bash
conda activate hybas
export CONDARC="$PWD/.condarc"
snakemake -s Hybas.smk --use-conda --cores 8
```

Use `-n` to inspect the DAG, `--printshellcmds` for command details, or `--keep-going` to finish independent isolates after a failure. The default profile enables `--rerun-incomplete` and displays logs for failed jobs.

The input layout, polishing route table, database preparation commands and output paths are documented in the [main README](../README.md). Final output is `<base_dir>/<isolate>/work/final/assembly.final.fasta`. Consult `metadata/polish_execution.json` for the actual executed stages rather than inferring them from `polish_mode`.

A report refresh does not require forcing assembly or Kraken. To refresh summaries explicitly:

```bash
snakemake -s Hybas.smk --use-conda --cores 8 --forcerun run_summary
```

MultiQC depends on those summaries and refreshes downstream. The optional UI launches the same Snakemake workflow; no additional launcher is needed.

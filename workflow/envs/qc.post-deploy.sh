#!/usr/bin/env bash
set -euo pipefail
# Run inside Snakemake's environment; never activate another environment here.
for tool in fastqc NanoPlot kraken2 extract_kraken_reads.py busco quast.py multiqc seqkit samtools pigz java; do
  location=$(command -v "$tool") || { echo "Missing executable: $tool" >&2; exit 1; }
  case "$location" in "$CONDA_PREFIX"/*) ;; *) echo "Executable outside managed environment: $tool" >&2; exit 1;; esac
done
java -version >/dev/null 2>&1
samtools --version >/dev/null

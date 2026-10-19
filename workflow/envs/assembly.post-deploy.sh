#!/usr/bin/env bash
set -euo pipefail
# Run inside Snakemake's environment; never activate another environment here.
for tool in unicycler racon minimap2 fastp filtlong pigz samtools spades.py; do
  location=$(command -v "$tool") || { echo "Missing executable: $tool" >&2; exit 1; }
  case "$location" in "$CONDA_PREFIX"/*) ;; *) echo "Executable outside managed environment: $tool" >&2; exit 1;; esac
done
samtools --version >/dev/null

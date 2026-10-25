#!/usr/bin/env bash
set -euo pipefail
# Run inside Snakemake's environment; never activate another environment here.
for tool in bwa samtools polypolish pigz; do
  location=$(command -v "$tool") || { echo "Missing executable: $tool" >&2; exit 1; }
  case "$location" in "$CONDA_PREFIX"/*) ;; *) echo "Executable outside managed environment: $tool" >&2; exit 1;; esac
done
samtools --version >/dev/null
polypolish --version >/dev/null

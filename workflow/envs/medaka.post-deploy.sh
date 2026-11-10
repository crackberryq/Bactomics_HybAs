#!/usr/bin/env bash
set -euo pipefail
# Run inside Snakemake's environment; never activate another environment here.
for tool in medaka_consensus medaka; do
  location=$(command -v "$tool") || { echo "Missing executable: $tool" >&2; exit 1; }
  case "$location" in "$CONDA_PREFIX"/*) ;; *) echo "Executable outside managed environment: $tool" >&2; exit 1;; esac
done
medaka --version >/dev/null

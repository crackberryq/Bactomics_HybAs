#!/usr/bin/env bash
set -euo pipefail
for tool in spades.py python pigz; do
  location=$(command -v "$tool") || { echo "Missing executable: $tool" >&2; exit 1; }
  case "$location" in "$CONDA_PREFIX"/*) ;; *) echo "Executable outside managed environment: $tool" >&2; exit 1;; esac
done
spades.py --version | grep -F '4.3.0'

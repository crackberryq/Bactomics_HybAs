#!/usr/bin/env bash
set -euo pipefail
for tool in flye python; do
  location=$(command -v "$tool")
  case "$location" in "$CONDA_PREFIX"/*) ;; *) echo "Executable outside managed environment: $tool" >&2; exit 1;; esac
done
flye --version | grep -F '2.9.6'

#!/usr/bin/env bash
set -euo pipefail
for tool in python fastp nxtrim; do
  location=$(command -v "$tool")
  case "$location" in "$CONDA_PREFIX"/*) ;; *) echo "Executable outside managed environment: $tool" >&2; exit 1;; esac
done
fastp --version
nxtrim 2>&1 | grep -F '0.4.3'

#!/usr/bin/env bash
set -euo pipefail
iso=$(realpath "${1:?Supply the completed isolate output directory}")
threads=${2:-4}
source_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace=${BACTOLOGIX_WORKSPACE:-"$HOME/.local/share/bactologix-dev"}
conda=${CONDA_EXE:-"$workspace/tools/miniforge/bin/conda"}
spec="$source_dir/workflow/envs/polishing_validation.yml"
digest=$(sha256sum "$spec" | cut -d' ' -f1)
prefix="$workspace/envs/polishing-validation/$digest"
mkdir -p "$(dirname "$prefix")"
exec 9>"$prefix.lock"
flock 9
if [[ ! -f "$prefix.ready" ]]; then
  "$conda" env create -p "$prefix" -f "$spec" -y
  touch "$prefix.ready"
fi
flock -u 9
cd "$iso/reports"
"$conda" run --no-capture-output -p "$prefix" python "$source_dir/workflow/scripts/polishing_validation.py" --isolate-dir "$iso" --threads "$threads" --refresh-report

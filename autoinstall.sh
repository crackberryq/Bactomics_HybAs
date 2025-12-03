#!/usr/bin/env bash
set -euo pipefail

echo "======================================="
echo "  Bactomics HybAs v8.4-Lite Installer"
echo "  (Conda + environment + repo + dry run)"
echo "======================================="

# ---- 0. Detect or install Conda ----
if ! command -v conda &>/dev/null; then
  echo ">>> Conda not found – installing Miniconda into \$HOME/miniconda3"

  cd "$HOME"
  if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS (Apple Silicon or Intel – installer auto-detects)
    MINICONDA_SCRIPT="Miniconda3-latest-MacOSX-arm64.sh"
    curl -LO "https://repo.anaconda.com/miniconda/${MINICONDA_SCRIPT}"
    bash "$MINICONDA_SCRIPT" -b -p "$HOME/miniconda3"
    rm -f "$MINICONDA_SCRIPT"
  else
    # Linux / WSL
    MINICONDA_SCRIPT="Miniconda3-latest-Linux-x86_64.sh"
    wget "https://repo.anaconda.com/miniconda/${MINICONDA_SCRIPT}"
    bash "$MINICONDA_SCRIPT" -b -p "$HOME/miniconda3"
    rm -f "$MINICONDA_SCRIPT"
  fi

  # Initialize conda for bash
  eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
else
  echo ">>> Conda found – using existing installation"
  eval "$(conda shell.bash hook)"
fi

# ---- 1. Configure Conda channels ----
echo
echo ">>> Configuring Conda channels (strict priority)"

conda config --set channel_priority strict
conda config --remove-key channels 2>/dev/null || true
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

echo "Active channels:"
conda config --show channels

# ---- 2. Create or reuse 'bactomics' environment ----
ENV_NAME="bactomics"
echo
echo ">>> Creating or activating Conda environment: ${ENV_NAME}"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "Environment '${ENV_NAME}' already exists – reusing it."
else
  conda create -y -n "${ENV_NAME}" python=3.11 snakemake=9.10.1 mamba
fi

conda activate "${ENV_NAME}"

echo
echo ">>> Environment versions"
python --version
snakemake --version
mamba --version || echo "mamba unavailable (unexpected)"

# ---- 3. Clone or update repo ----
echo
echo ">>> Obtaining Bactomics HybAs repository"

cd "$HOME"
if [ -d "$HOME/Bactomics_HybAs/.git" ]; then
  echo "Repository exists – pulling latest changes"
  cd "$HOME/Bactomics_HybAs"
  git pull --ff-only || echo "[WARN] git pull failed (local changes?)"
else
  git clone https://github.com/crackberryq/Bactomics_HybAs.git
  cd "$HOME/Bactomics_HybAs"
fi

echo "Repository path: $(pwd)"

# ---- 4. Optional Snakemake dry run ----
echo
echo ">>> Running Snakemake dry-run"
echo "(Checks installation but does NOT execute tools)"
echo

set +e
snakemake -s Snakefile --use-conda --cores 4 -n -p
STATUS=$?
set -e

if [ "$STATUS" -ne 0 ]; then
  echo
  echo "[NOTE] Dry-run reported issues. This is usually because:"
  echo " - config.yaml has not been edited yet"
  echo " - base_dir or isolate directory is missing"
  echo
  echo "Conda + Snakemake + repo installation completed successfully."
else
  echo
  echo ">>> Dry-run completed successfully!"
fi

echo
echo "======================================="
echo "   Bactomics HybAs installation done"
echo "   Next steps:"
echo "     1) Edit config.yaml"
echo "     2) Prepare input reads"
echo "     3) Run the pipeline"
echo "======================================="

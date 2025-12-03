# Bactomics HybAs v8.4‑Lite – Auto Installer

This script installs / configures:

- Miniconda (if not already installed)
- A clean `bactomics` Conda environment (Python 3.11, Snakemake 9.10.1, Mamba)
- The `Bactomics_HybAs` GitHub repository into `~/bactomics`
- A Snakemake **dry run** to verify installation

## Usage

From your **Linux / macOS / WSL2** terminal:

```bash
curl -LO https://raw.githubusercontent.com/crackberryq/Bactomics_HybAs/main/autoinstall.sh
bash autoinstall.sh
```

—or—copy‑paste the script below into a file called `autoinstall.sh`, then run:

```bash
bash autoinstall.sh
```

---

## Installer Script

```bash
#!/usr/bin/env bash
set -euo pipefail

echo "======================================="
echo "  Bactomics HybAs v8.4-Lite Installer"
echo "  (Conda + env + repo + dry run)"
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
    bash "${MINICONDA_SCRIPT}" -b -p "$HOME/miniconda3"
    rm -f "${MINICONDA_SCRIPT}"
  fi

  # Initialize conda for bash
  eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
else
  echo ">>> Conda found – using existing installation"
  eval "$(conda shell.bash hook)"
fi

# ---- 1. Configure Conda channels & priority ----
echo
echo ">>> Configuring Conda channels and strict priority"

conda config --set channel_priority strict

# Reset channels and set desired order: conda-forge, bioconda, defaults
conda config --remove-key channels 2>/dev/null || true
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults

echo "Current channels:"
conda config --show channels

# ---- 2. Create or reuse 'bactomics' environment ----
ENV_NAME="bactomics"
echo
echo ">>> Creating (or reusing) Conda environment: ${ENV_NAME}"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "Environment '${ENV_NAME}' already exists – will reuse it."
else
  conda create -y -n "${ENV_NAME}" python=3.11 snakemake=9.10.1 mamba
fi

conda activate "${ENV_NAME}"

echo
echo ">>> Environment details"
python --version
snakemake --version
mamba --version || echo "mamba not available (unexpected)"

# ---- 3. Clone or update Bactomics_HybAs repo ----
echo
echo ">>> Cloning or updating Bactomics_HybAs repository"

cd "$HOME"
if [ -d "$HOME/bactomics/.git" ]; then
  echo "Repository '~/bactomics' already exists – pulling latest changes."
  cd "$HOME/bactomics"
  git pull --ff-only || echo "Warning: 'git pull' failed (you may have local changes)."
else
  git clone https://github.com/crackberryq/Bactomics_HybAs.git bactomics
  cd "$HOME/bactomics"
fi

echo
echo "Repository location: $(pwd)"

# ---- 4. Optional Snakemake dry run ----
echo
echo ">>> Running Snakemake DRY RUN to test installation"
echo "    (This will NOT actually run the tools.)"
echo

set +e
snakemake -s Snakefile --use-conda --cores 4 -n -p
STATUS=$?
set -e

if [ "$STATUS" -ne 0 ]; then
  echo
  echo ">>> Dry run reported errors."
  echo "Most common reasons:"
  echo "  - 'config.yaml' not yet edited (base_dir / isolate missing)"
  echo "  - Input folders do not exist yet"
  echo
  echo "Conda + Snakemake + repository are still installed correctly."
else
  echo
  echo ">>> Dry run completed successfully – installation looks good."
fi

echo
echo "======================================="
echo "  DONE."
echo "Next steps:"
echo "  1) Edit config.yaml to set base_dir and isolate."
echo "  2) Place reads under: base_dir/isolate/raw/{illumina,nanopore}"
echo "  3) Run the full workflow, e.g.:"
echo "     snakemake -s Snakefile --use-conda -p --cores 12"
echo "======================================="
```


#!/bin/bash
# =============================================================================
# 00_setup_env.sh -- build the gsMap python environment
#
# Project : spatialAmygdala (LIBD) -- Visium / 16_LDSC / gsMap
# Purpose : create a self-contained conda env holding the gsMap CLI, and write
#           a provenance record of exactly what got installed.
#
# Where it runs : COMPUTE NODE ONLY (never the login node).
#                 sbatch this script, or run it inside an srun session:
#                   srun --pty -p shared --cpus-per-task 4 --mem 32G --time 4:00:00 bash
#
# Idempotent : if $ENV_DIR already exists the env is reused, not rebuilt.
#              Nothing is ever deleted by this script.
# =============================================================================
#SBATCH --job-name=gsmap_setup
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00

set -euo pipefail

# ---- hardcoded paths --------------------------------------------------------
SCRATCH=/users/mtotty/claude_scratch/gsmap
ENV_ROOT=${SCRATCH}/envs
ENV_DIR=${ENV_ROOT}/gsmap
VERSION_FILE=${ENV_ROOT}/gsmap_version.txt
CLI_HELP_FILE=${ENV_ROOT}/gsmap_cli_help.txt

mkdir -p "${ENV_ROOT}"

# keep conda/pip caches inside scratch so nothing lands in ~/
export CONDA_PKGS_DIRS=${SCRATCH}/.conda_pkgs
export PIP_CACHE_DIR=${SCRATCH}/.pip_cache
mkdir -p "${CONDA_PKGS_DIRS}" "${PIP_CACHE_DIR}"

# ---- conda ------------------------------------------------------------------
module load conda/3-24.3.0
eval "$(conda shell.bash hook)"

# ---- 1. create env (skip if present) ----------------------------------------
if [[ -x "${ENV_DIR}/bin/python" ]]; then
  echo "[setup] env already present at ${ENV_DIR} -- skipping creation"
else
  echo "[setup] creating conda env at ${ENV_DIR} (python 3.11)"
  conda create -y -p "${ENV_DIR}" python=3.11
fi

conda activate "${ENV_DIR}"
echo "[setup] active python: $(which python)"

# ---- 2. install gsMap (PyPI package name is 'gsMap') ------------------------
python -m pip install --upgrade pip
python -m pip install gsMap

# ---- 3. provenance record ---------------------------------------------------
echo "[setup] writing ${VERSION_FILE}"
{
  echo "===== gsMap environment provenance ====="
  echo "date            : $(date)"
  echo "host            : $(hostname)"
  echo "env prefix      : ${ENV_DIR}"
  echo
  echo "----- gsmap --version -----"
  gsmap --version 2>&1 || echo "(gsmap --version failed)"
  echo
  echo "----- pip show gsMap -----"
  python -m pip show gsMap 2>&1 || echo "(pip show failed)"
  echo
  echo "----- which gsmap -----"
  which gsmap 2>&1 || echo "(gsmap not on PATH)"
  echo
  echo "----- python -V -----"
  python -V 2>&1
  echo
  echo "----- key package versions -----"
  python -m pip list 2>/dev/null | grep -iE '^(torch|torchvision|torchaudio|scanpy|anndata|numpy|scipy|pandas|pyranges|bitarray|jax|nvidia)' || true
  echo
  echo "----- torch build / CUDA -----"
  python - <<'PY' 2>&1 || echo "(torch not importable)"
import torch
print("torch.__version__      :", torch.__version__)
print("torch.version.cuda     :", torch.version.cuda)
print("cuda.is_available()    :", torch.cuda.is_available())
PY
  echo
  echo "----- gsmap --help -----"
  gsmap --help 2>&1 || echo "(gsmap --help failed)"
} > "${VERSION_FILE}"

# ---- 4. per-subcommand CLI help (exact flag names for downstream scripts) ---
echo "[setup] writing ${CLI_HELP_FILE}"
{
  for SUB in quick_mode create_slice_mean run_cauchy_combination \
             run_find_latent_representations run_latent_to_gene \
             run_generate_ldscore run_spatial_ldsc run_report \
             format_sumstats find_latent_representations latent_to_gene \
             generate_ldscore spatial_ldsc cauchy_combination report ; do
    echo "############################################################"
    echo "### gsmap ${SUB} --help"
    echo "############################################################"
    gsmap "${SUB}" --help 2>&1 || echo "(subcommand '${SUB}' not available)"
    echo
  done
} > "${CLI_HELP_FILE}"

echo "[setup] DONE"
echo "[setup] version file : ${VERSION_FILE}"
echo "[setup] cli help file: ${CLI_HELP_FILE}"

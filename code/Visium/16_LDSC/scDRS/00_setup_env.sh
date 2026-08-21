#!/bin/bash
#SBATCH --job-name=scdrs_env
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=2:00:00
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/00_setup_env.%j.log
# -----------------------------------------------------------------------------
# 00_setup_env.sh -- build the scDRS conda env in scratch. Run ONCE.
#
#   sbatch 00_setup_env.sh
#
# Env creation is CPU/RAM heavy -- never run it on the login node.
# Everything lands in /users/mtotty/claude_scratch/scdrs/envs/scdrs.
# Nothing outside scratch is written.
# -----------------------------------------------------------------------------
set -euo pipefail
# NOTE: under sbatch, SLURM copies this script to a spool directory, so "$0"
# does NOT resolve to this file's real location -- $(dirname $(readlink -f $0))
# yields /var/spool/slurm/... and config.sh is not found there. Hardcode the
# code directory instead. Override with CODE_DIR=... if you relocate the tree.
CODE_DIR=${CODE_DIR:-/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS}
source "$CODE_DIR/config.sh"

echo "**** Job starts ****"; date; hostname

mkdir -p "$SCRATCH" "$SCRATCH/envs" "$CODE/logs"

# keep conda/pip caches out of ~/ (home quota is small)
export CONDA_PKGS_DIRS=$SCRATCH/.conda_pkgs
export PIP_CACHE_DIR=$SCRATCH/.pip_cache
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR"

module load $CONDA_MODULE
eval "$(conda shell.bash hook)"

if [ -x "$SCDRS" ]; then
    echo "env already exists at $ENVDIR -- skipping create"
else
    conda create -y -p "$ENVDIR" python=3.11
fi

conda activate "$ENVDIR"

# scDRS 1.0.2 is the current PyPI release. Its deps (scanpy, anndata,
# scikit-misc, statsmodels, fire) come along automatically.
pip install --no-cache-dir "scdrs==1.0.2"
pip install --no-cache-dir matplotlib seaborn pyarrow

echo "---- versions ----"
"$PY" - <<'PY'
import scdrs, scanpy as sc, anndata, numpy, pandas, scipy
print("scdrs   ", scdrs.__version__)
print("scanpy  ", sc.__version__)
print("anndata ", anndata.__version__)
print("numpy   ", numpy.__version__)
print("pandas  ", pandas.__version__)
print("scipy   ", scipy.__version__)
PY

echo "---- CLI smoke test ----"
"$SCDRS" --help | head -30

echo
echo "Activate later with:"
echo "  module load $CONDA_MODULE"
echo "  eval \"\$(conda shell.bash hook)\""
echo "  conda activate $ENVDIR"

echo "**** Job ends ****"; date

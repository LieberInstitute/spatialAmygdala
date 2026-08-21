#!/bin/bash
#SBATCH --job-name=scdrs_controls
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/09_controls.%j.log
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/09_controls.%j.log
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#
# Control analyses for the scDRS results paragraph.
#
#   09_controls_compute.py  -> ctrl_*.tsv   (statistics; the slow part)
#   10_controls_figures.py  -> fig_ctrl_*.png
#
# Re-running only the figures after editing style:
#     sbatch --export=ALL,FIGURES_ONLY=1 09_run_controls.sh
#
# Memory: reads 42 score files one at a time (~300 MB each, one column
# retained) plus the h5ad obs table read in backed mode -- peak stays well
# under 32 G. Runtime is dominated by the 42 sequential score reads.
#
# Writes ONLY into $WORKDIR/controls/. Nothing is deleted or overwritten
# elsewhere.

set -euo pipefail

CODE_DIR="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS"
ENVDIR="/users/mtotty/claude_scratch/scdrs/envs/scdrs"
FIGURES_ONLY="${FIGURES_ONLY:-0}"

echo "=== scDRS controls ==="
echo "host      : $(hostname)"
echo "started   : $(date)"
echo "code dir  : $CODE_DIR"
echo "figs only : $FIGURES_ONLY"

module load conda/3-24.3.0
eval "$(conda shell.bash hook)"
conda activate "$ENVDIR"
echo "python    : $(which python)"

cd "$CODE_DIR"

if [ "$FIGURES_ONLY" != "1" ]; then
    echo
    echo "--- 09_controls_compute.py ---"
    python 09_controls_compute.py
else
    echo "(skipping compute step; FIGURES_ONLY=1)"
fi

echo
echo "--- 10_controls_figures.py ---"
python 10_controls_figures.py

echo
echo "=== outputs ==="
ls -la /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/scDRS/controls/

echo
echo "finished  : $(date)"

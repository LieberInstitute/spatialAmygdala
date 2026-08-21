#!/bin/bash
#SBATCH --job-name=scdrs_figs
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/07_figures.%j.log
# -----------------------------------------------------------------------------
# 07_run_figures.sh -- render figures from whatever scores exist.
#
#   sbatch 07_run_figures.sh                       # everything under score/
#   SUB=pilot sbatch 07_run_figures.sh             # just the pilot outputs
# -----------------------------------------------------------------------------
set -euo pipefail
# NOTE: under sbatch, SLURM copies this script to a spool directory, so "$0"
# does NOT resolve to this file's real location -- $(dirname $(readlink -f $0))
# yields /var/spool/slurm/... and config.sh is not found there. Hardcode the
# code directory instead. Override with CODE_DIR=... if you relocate the tree.
CODE_DIR=${CODE_DIR:-/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS}
source "$CODE_DIR/config.sh"

echo "**** Job starts ****"; date; hostname
activate_env

SUB=${SUB:-}
if [ -n "$SUB" ]; then
    SCORE_SUB="score/$SUB"; DOWN_SUB="downstream/$SUB"; OUT="$FIG_DIR/$SUB"
else
    SCORE_SUB="score"; DOWN_SUB="downstream"; OUT="$FIG_DIR"
fi
mkdir -p "$OUT" "$CODE/logs"

"$PY" "$CODE/07_figures.py" \
    --workdir "$WORKDIR" \
    --h5ad "$H5AD_SCDRS" \
    --annot "$ANNOT_COL" \
    --ldsc "$PROJ/code/Visium/16_LDSC/LDSC/ldsc_results.csv" \
    --outdir "$OUT" \
    --score-subdir "$SCORE_SUB" \
    --downstream-subdir "$DOWN_SUB"

ls -lh "$OUT"
echo "**** Job ends ****"; date

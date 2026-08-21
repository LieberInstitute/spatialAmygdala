#!/bin/bash
#SBATCH --job-name=scdrs_h5ad
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/04_prep_h5ad.%j.log
# -----------------------------------------------------------------------------
# 04_prep_h5ad.sh -- concatenate the seven gsMap section h5ad files into the
# single scDRS input, and write the covariate table.
#
#   sbatch 04_prep_h5ad.sh
#
# Reads (read-only):  ../gsMap/ST/*.h5ad  -- 224,021 spots x 19,389 genes total
# Writes:             $H5AD_SCDRS, $COV_FILE
#
# RESOURCES: 4 CPU / 48 G. The seven section matrices total ~600 M non-zeros;
# held as float64 CSR that is ~7 GB, and concat transiently holds both the
# inputs and the result. 48 G covers that with roughly 3x headroom. Single-
# threaded work apart from gzip, so the 4 cores are for the compressed write.
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
mkdir -p "$(dirname "$H5AD_SCDRS")" "$CODE/logs"

if [ -f "$H5AD_SCDRS" ]; then
    echo "$H5AD_SCDRS already exists -- remove it yourself if you want a rebuild"
    echo "**** Job ends ****"; date; exit 0
fi

"$PY" "$CODE/04a_prep_h5ad.py" \
    "$ST_DIR" Br6471 Br6660 Br6423 Br2743 Br8325 Br9192 Br9280 \
    --out "$H5AD_SCDRS" \
    --cov "$COV_FILE" \
    --annot "$ANNOT_COL" \
    | tee "$CODE/04_h5ad_summary.txt"

ls -lh "$H5AD_SCDRS" "$COV_FILE"
echo "**** Job ends ****"; date

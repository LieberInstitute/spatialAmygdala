#!/bin/bash
#SBATCH --job-name=scdrs_contrib
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/06e_domain_gene_contributions.%j.log

set -euo pipefail

CODE_DIR=${CODE_DIR:-/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS}
source "$CODE_DIR/config.sh"

echo "**** Job starts ****"
date
hostname

activate_env

mkdir -p "$FIG_DIR" "$CODE_DIR/logs"

GS=$GS_DIR/magma_top${GS_NMAX}.gs

/usr/bin/time -v "$PY" "$CODE_DIR/06d_domain_gene_contributions.py" \
    "$H5AD_SCDRS" \
    "$COV_FILE" \
    "$GS" \
    "$FIG_DIR"

echo
echo "---- outputs ----"

ls -lh \
    "$FIG_DIR/SCZ_LA_gene_contributions.tsv" \
    "$FIG_DIR/MDD_BM_gene_contributions.tsv" \
    "$FIG_DIR/BIP_2024_CLA_gene_contributions.tsv"

echo
echo "**** Job ends ****"
date
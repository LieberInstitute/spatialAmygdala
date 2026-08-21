#!/bin/bash
#SBATCH --job-name=scdrs_psych_donor
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/06g_donor_psych_enrichment.%j.log

set -euo pipefail

CODE_DIR=${CODE_DIR:-/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS}
source "$CODE_DIR/config.sh"

echo "**** Job starts ****"
date
hostname

activate_env

mkdir -p "$FIG_DIR" "$CODE_DIR/logs"

/usr/bin/time -v "$PY" "$CODE_DIR/06f_donor_psych_enrichment.py" \
    "$H5AD_SCDRS" \
    "$SCORE_DIR" \
    "$FIG_DIR"

echo
echo "---- outputs ----"

ls -lh \
    "$FIG_DIR/psych_enrichment_donor_trait_domain.tsv" \
    "$FIG_DIR/psych_enrichment_donor_domain.tsv" \
    "$FIG_DIR/psych_enrichment_domain_summary.tsv"

echo
echo "**** Job ends ****"
date
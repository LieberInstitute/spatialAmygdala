#!/bin/bash
#SBATCH --job-name=scdrs_down
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=18:00:00
#SBATCH --array=1-42%20
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/06c_down.%A_%a.log
# -----------------------------------------------------------------------------
# 06c_downstream_array.sh -- group analysis per trait, as an array.
#
#   sbatch 06c_downstream_array.sh   # after 06b has produced .full_score.gz
#
# RUNTIME. The pilot suggested ~15 min per trait; the real spread on the full
# run (job 34733151) was 20 min to >7 h. Three traits (EduYears, Smoking,
# Stroke_2022_Any) hit a 6 h wall and were killed with nothing written, so the
# limit is now 18 h. Cost tracks the number of spots passing the FDR filter,
# not the trait's rank in the .gs file -- traits with broad enrichment have
# more spots to permute in the heterogeneity test.
#
# Memory: this step loads the same 224k x 19k matrix as 06b plus one trait's
# full_score, so it sits in the same ~14 GB class. 24 G matches 06b.
#
# REQUIRES the NumPy 2.x patch to scdrs/method.py -- see README. Without it
# this step dies with "np.float_ was removed in the NumPy 2.0 release".
#
# IDEMPOTENT: skips a trait whose group file already exists. Never deletes.
# -----------------------------------------------------------------------------
set -euo pipefail

CODE_DIR=${CODE_DIR:-/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS}
source "$CODE_DIR/config.sh"

echo "**** Job starts ****"; date; hostname
activate_env
mkdir -p "$DOWNSTREAM_DIR" "$CODE/logs"

GS=${GS_CHUNK:-$GS_DIR/magma_top${GS_NMAX}.gs}
LINE=$(( ${SLURM_ARRAY_TASK_ID:-1} + 1 ))
TRAIT=$(sed -n "${LINE}p" "$GS" | cut -f1)

if [ -z "${TRAIT:-}" ]; then
    echo "no trait at .gs line $LINE -- array index out of range"; exit 0
fi
echo "task ${SLURM_ARRAY_TASK_ID:-1} -> trait '$TRAIT'"

OUT=$DOWNSTREAM_DIR/$TRAIT.scdrs_group.$ANNOT_COL
if [ -s "$OUT" ]; then
    echo "$TRAIT group analysis already present -- nothing to do"; ls -lh "$OUT"
    echo "**** Job ends ****"; date; exit 0
fi

if [ ! -s "$SCORE_DIR/$TRAIT.full_score.gz" ]; then
    echo "ERROR: $SCORE_DIR/$TRAIT.full_score.gz missing -- run 06b first"; exit 1
fi

# --score-file takes a pattern with @ standing for the trait name; giving the
# literal trait restricts this task to exactly one trait.
/usr/bin/time -v "$SCDRS" perform-downstream \
    --h5ad-file "$H5AD_SCDRS" \
    --score-file "$SCORE_DIR/$TRAIT.full_score.gz" \
    --out-folder "$DOWNSTREAM_DIR" \
    --group-analysis "$ANNOT_COL" \
    --flag-filter-data True \
    --flag-raw-count True

ls -lh "$DOWNSTREAM_DIR/$TRAIT".*
echo "**** Job ends ****"; date

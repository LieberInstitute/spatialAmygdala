#!/bin/bash
#SBATCH --job-name=scdrs_pilot
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/05_score_pilot.%j.log
# -----------------------------------------------------------------------------
# 05_score_pilot.sh -- PILOT. Score three traits (SCZ, MDD, Height) and run the
# domain-level group analysis, so we learn the real runtime and memory profile
# before committing to all 40.
#
#   sbatch 05_score_pilot.sh
#
# SCZ and MDD are the positive controls (both are significant in the classic
# s-LDSC results and in the gsMap pilot); Height is the negative control and
# should come out flat across every domain. If Height lights up, the covariate
# model or the gene sets are wrong and the sweep should not be launched.
#
# WHY THESE FLAGS
#   --flag-filter-data True    drop genes/spots below scDRS's own thresholds
#   --flag-raw-count True      X holds raw counts, so scDRS does its own
#                              normalisation + log1p internally
#   --n-ctrl 1000              control gene sets; 1000 is the paper's default
#                              and gives p-values down to ~1e-3 empirically,
#                              with the MC-adjusted tail extending further
#   --cov-file                 regress out depth + donor (see 04a_prep_h5ad.py)
#   --flag-return-ctrl-raw-score False / --flag-return-ctrl-norm-score True
#                              keeps .full_score.gz to a workable size while
#                              retaining what group-analysis needs
#
# RESOURCES: 8 CPU / 48 G. Verified against the scdrs 1.0.2 source: when .X is
# sparse AND a --cov-file is given, preprocess uses IMPLICIT covariate
# correction (scdrs/pp.py) -- it stores the covariate betas and applies them
# chunk-wise instead of materialising the corrected dense matrix. The 224k x
# 19k matrix therefore stays sparse (~600 M non-zeros, ~7 GB as float32 CSR).
# A dense copy would have been ~17 GB, which is why the first draft of this
# script asked for 192 G. 48 G leaves room for the control-score arrays.
# Watch the /usr/bin/time -v "Maximum resident set size" line below and
# right-size 06_score_full.sh from it.
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
mkdir -p "$SCORE_DIR/pilot" "$DOWNSTREAM_DIR/pilot" "$CODE/logs"

GS=$GS_DIR/magma_top${GS_NMAX}.gs
PILOT_GS=$GS_DIR/pilot.gs

# subset the .gs to the three pilot traits (header line + the three rows)
awk 'NR==1 || $1=="SCZ" || $1=="MDD" || $1=="Height"' "$GS" > "$PILOT_GS"
echo "pilot gene sets:"; cut -f1 "$PILOT_GS"

echo
echo "[1/2] scdrs compute-score"
/usr/bin/time -v "$SCDRS" compute-score \
    --h5ad-file "$H5AD_SCDRS" \
    --h5ad-species human \
    --gs-file "$PILOT_GS" \
    --gs-species human \
    --cov-file "$COV_FILE" \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl "$N_CTRL" \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    --out-folder "$SCORE_DIR/pilot"

ls -lh "$SCORE_DIR/pilot"

echo
echo "[2/2] scdrs perform-downstream (group analysis by $ANNOT_COL)"
"$SCDRS" perform-downstream \
    --h5ad-file "$H5AD_SCDRS" \
    --score-file "$SCORE_DIR/pilot/@.full_score.gz" \
    --out-folder "$DOWNSTREAM_DIR/pilot" \
    --group-analysis "$ANNOT_COL" \
    --flag-filter-data True \
    --flag-raw-count True

echo "---- group analysis results ----"
for f in "$DOWNSTREAM_DIR"/pilot/*.scdrs_group.*; do
    echo "== $f"; cat "$f"
done

echo "**** Job ends ****"; date

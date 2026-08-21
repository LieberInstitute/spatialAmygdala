#!/bin/bash
#SBATCH --job-name=cond_cauchy_x
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --exclude=compute-058,compute-175,compute-053
# -----------------------------------------------------------------------------
# 13_cond_cauchy_across.sh -- ONE trait: domain-level p-values pooled ACROSS all
# 7 capture areas, for the CONDITIONED arm. Mirrors 05_cauchy_across_samples.sh
# (which does this for arm 0 / unconditioned) with three deliberate changes:
#   * one trait per job instead of a 40-trait serial loop, so the whole stage
#     fits the 1-day walltime cap and a single bad trait cannot strand the rest;
#   * --workdir points at gsMap_cond_${ARM};
#   * TRAIT arrives via --export.
#
# As in 05, pooling across samples is only meaningful because every sample was
# ranked against the SHARED slice mean, and Cauchy/ACAT is used because spot-level
# p-values within a domain are strongly dependent (neighbouring spots share
# LD-score signal), which Fisher's method would not tolerate.
#
#   sbatch --export=ALL,ARM=functional,TRAIT=SCZ 13_cond_cauchy_across.sh
#
# OUTPUT: $WORKDIR/cauchy_across_samples/<trait>_cauchy.csv.gz
# -----------------------------------------------------------------------------
set -uo pipefail
ARM="${ARM:?}"; TRAIT="${TRAIT:?}"
CODE=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
WORKDIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}
ANNOT=BS_k16_Semisupervised_wAI
OUTDIR=$WORKDIR/cauchy_across_samples
export PATH="/users/mtotty/claude_scratch/gsmap/envs/gsmap/bin:$PATH"
mkdir -p "$OUTDIR"

echo "**** cauchy start ****"; date
echo "ARM=$ARM TRAIT=$TRAIT host=$(hostname) job=${SLURM_JOB_ID:-NA}"

SAMPLE_LIST=$(tr '\n' ' ' < "$CODE/samples.txt")
# every one of the 7 per-spot inputs must exist, or the pooled p-value would be
# silently computed over a subset of sections.
MISS=0
for S in $SAMPLE_LIST; do
  F="$WORKDIR/$S/spatial_ldsc/${S}_${TRAIT}.csv.gz"
  [ -s "$F" ] || { echo "MISSING INPUT: $F"; MISS=$((MISS+1)); }
done
[ "$MISS" -gt 0 ] && { echo "FATAL: $MISS/7 per-spot inputs absent for $TRAIT"; exit 7; }

OUT=$OUTDIR/${TRAIT}_cauchy.csv.gz
gsmap run_cauchy_combination \
  --workdir "$WORKDIR" \
  --trait_name "$TRAIT" \
  --annotation "$ANNOT" \
  --sample_name_list $SAMPLE_LIST \
  --output_file "$OUT"
RC=$?
echo "gsmap rc=$RC (non-fatal; success judged by output file)"
if [ -s "$OUT" ]; then echo "OK $TRAIT -> $OUT"; ls -l "$OUT"
else echo "FAIL: no cross-sample Cauchy output for $TRAIT (rc=$RC)"; exit 6; fi
echo "**** cauchy end ****"; date

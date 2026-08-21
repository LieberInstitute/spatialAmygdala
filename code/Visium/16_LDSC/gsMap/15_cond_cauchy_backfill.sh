#!/bin/bash
#SBATCH --job-name=cond_cauchy_bf
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --exclude=compute-058,compute-175,compute-053
# -----------------------------------------------------------------------------
# 15_cond_cauchy_backfill.sh <ARM> -- per-sample Cauchy for units whose per-spot
# regression already completed but whose Cauchy output is absent.
#
# WHY THIS EXISTS: run_cauchy_combination reads
#   {workdir}/<S>/find_latent_representations/<S>_add_latent.h5ad
# which was never carried into the conditional workdirs, so the Cauchy step died
# with FileNotFoundError while the expensive 2.5-3 h regression had succeeded.
# The h5ad is now symlinked from the arm-0 workdir. This backfills the gap
# WITHOUT recomputing any regression (~3 s per unit in the pilot).
#
# Re-runnable: skips any unit whose Cauchy output already exists.
# -----------------------------------------------------------------------------
set -uo pipefail
ARM="${ARM:-functional}"
CODE=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
WORKDIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}
export PATH="/users/mtotty/claude_scratch/gsmap/envs/gsmap/bin:$PATH"
echo "**** backfill start ****"; date; echo "host=$(hostname) ARM=$ARM"

nok=0; nfail=0; nskip=0
while read -r S; do
  [ -z "$S" ] && continue
  H="$WORKDIR/$S/find_latent_representations/${S}_add_latent.h5ad"
  [ -s "$H" ] || { echo "FATAL: latent h5ad missing/unreadable for $S: $H"; exit 8; }
  for F in "$WORKDIR/$S/spatial_ldsc/${S}_"*.csv.gz; do
    [ -e "$F" ] || continue
    B=$(basename "$F" .csv.gz); T=${B#${S}_}
    OUT="$WORKDIR/$S/cauchy_combination/${S}_${T}.Cauchy.csv.gz"
    if [ -s "$OUT" ]; then nskip=$((nskip+1)); continue; fi
    gsmap run_cauchy_combination --workdir "$WORKDIR" --sample_name "$S" \
      --trait_name "$T" --annotation BS_k16_Semisupervised_wAI >/dev/null 2>&1
    RC=$?
    if [ -s "$OUT" ]; then echo "OK   $S $T"; nok=$((nok+1))
    else echo "FAIL $S $T (rc=$RC)"; nfail=$((nfail+1)); fi
  done
done < "$CODE/samples.txt"
echo "----"; echo "backfilled=$nok failed=$nfail already-present=$nskip"
echo "**** backfill end ****"; date

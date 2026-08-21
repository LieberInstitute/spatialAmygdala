#!/bin/bash
#SBATCH --job-name=cond_topup
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --exclude=compute-058,compute-175,compute-053
# -----------------------------------------------------------------------------
# 16_cond_topup.sh -- keep the unit queue topped up and fire cross-sample Cauchy
# for traits that reach 7/7, then re-submit itself.
#
# Runs on a compute node (never a login node) and only calls squeue/sbatch.
#
# TIMEOUT BUG, fixed: the first version requested 2 h and slept 7 x 15 min = 1h45m
# plus per-pass overhead. It hit TIMEOUT at 02:00:21 BEFORE reaching its own
# re-submit line, so the self-perpetuating chain died silently and the sweep sat
# untended. Now: 4 h walltime, 6 passes (1h30m of sleeping), and the successor is
# submitted BEFORE the last sleep so a timeout can never break the chain again.
# -----------------------------------------------------------------------------
set -uo pipefail
ARM="${ARM:-functional}"; TARGET="${TARGET:-60}"
CODE=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
WORKDIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}
echo "topup start $(date) host=$(hostname) ARM=$ARM TARGET=$TARGET"

count_missing() {
  local n=0
  while read -r T; do while read -r S; do
    [ -s "$WORKDIR/$S/spatial_ldsc/${S}_${T}.csv.gz" ] || n=$((n+1))
  done < "$CODE/samples.txt"; done < "$CODE/traits.txt"
  echo "$n"
}

for i in $(seq 1 6); do
  NPEND=$(count_missing)
  echo "--- pass $i $(date) missing=$NPEND ---"
  if [ "$NPEND" -eq 0 ]; then
    bash "$CODE/17_cond_cauchy_sweep.sh" "$ARM" | tail -1
    echo "all units complete; chain ends after cauchy sweep"
    exit 0
  fi
  bash "$CODE/11_cond_submit_wave.sh" "$ARM" "$TARGET" | tail -2
  # fire cross-sample Cauchy incrementally as each trait reaches 7/7, so a
  # prerequisite failure surfaces immediately instead of after all 280 units
  bash "$CODE/17_cond_cauchy_sweep.sh" "$ARM" | tail -1

  # Submit the successor BEFORE the final sleep: if this job is killed or times
  # out during the sleep, the chain still survives.
  if [ "$i" -eq 6 ]; then
    NRUN=$(squeue -u "$USER" -h -o "%j" | grep -c "^cond_topup$")
    if [ "$NRUN" -gt 1 ]; then
      echo "successor already queued (n=$NRUN); not stacking another"
    else
      sbatch --export=ALL,ARM="$ARM",TARGET="$TARGET" \
        --output=$CODE/logs/cond_topup_%j.out --error=$CODE/logs/cond_topup_%j.out \
        "$CODE/16_cond_topup.sh"
      echo "re-queued self $(date)"
    fi
  fi
  sleep 900
done
echo "topup end $(date)"

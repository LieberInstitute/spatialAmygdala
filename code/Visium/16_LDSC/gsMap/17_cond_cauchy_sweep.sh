#!/bin/bash
# 17_cond_cauchy_sweep.sh <ARM> -- submit cross-sample Cauchy for every trait that
# has reached 7/7 per-spot outputs and does not yet have a pooled result.
#
# Run INCREMENTALLY (not batched at the end) so a broken prerequisite surfaces on
# the first eligible trait rather than after all 280 units finish. Idempotent:
# skips traits already done and traits already queued (squeue name lock).
#
# --time=01:00:00 deliberately: the step takes ~25 s, and a 1-day request sat
# behind the sweep's own units at equal priority instead of backfilling.
set -uo pipefail
ARM="${1:-functional}"
CODE=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
WORKDIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}
OUTDIR=$WORKDIR/cauchy_across_samples
mkdir -p "$OUTDIR"
INFLIGHT=$(squeue -u "$USER" -h -o "%j" | tr '\n' ' ')
nsub=0; nready=0; ndone=0
while read -r T; do
  [ -z "$T" ] && continue
  if [ -s "$OUTDIR/${T}_cauchy.csv.gz" ]; then ndone=$((ndone+1)); continue; fi
  n=0
  while read -r S; do
    [ -s "$WORKDIR/$S/spatial_ldsc/${S}_${T}.csv.gz" ] && n=$((n+1))
  done < "$CODE/samples.txt"
  [ "$n" -eq 7 ] || continue
  nready=$((nready+1))
  JN="cx_${ARM}_${T}"
  case " $INFLIGHT " in *" $JN "*) continue ;; esac
  JID=$(sbatch --parsable --job-name="$JN" --time=01:00:00 --cpus-per-task=2 --mem=48G \
    --export=ALL,ARM="$ARM",TRAIT="$T" \
    --output=$CODE/logs/${JN}_%j.out --error=$CODE/logs/${JN}_%j.out \
    "$CODE/13_cond_cauchy_across.sh")
  echo "submitted $JID $JN"; nsub=$((nsub+1))
done < "$CODE/traits.txt"
echo "cauchy: done=$ndone eligible-now=$nready submitted=$nsub"

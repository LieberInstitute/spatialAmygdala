#!/bin/bash
# 12_cond_monitor.sh <ARM> -- progress + stall detection for the per-unit sweep.
#
# sacct TotalCPU reads 00:00:00 for RUNNING jobs on this cluster, so it is useless
# for stall detection. We use log mtime age and the tqdm "Chunk-N/Total-chunk-37"
# progress line instead (tqdm writes carriage returns, hence tr '\r' '\n').
set -uo pipefail
ARM="${1:-functional}"
CODE=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
WORKDIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}

NOUT=$(cat $CODE/samples.txt | while read -r S; do ls -1 "$WORKDIR/$S/spatial_ldsc/"*.csv.gz 2>/dev/null; done | wc -l)
echo "OUTPUTS: $NOUT / 280"
echo "QUEUE: $(squeue -u $USER -h | wc -l) total, running=$(squeue -u $USER -h -t R | wc -l), pending=$(squeue -u $USER -h -t PD | wc -l)"

echo "--- RUNNING UNITS (chunk progress / log age) ---"
squeue -u "$USER" -h -t R -o "%i %j %M %N" | while read -r JID JN ELAP NODE; do
  LOG=$(ls -t $CODE/logs/${JN}_${JID}.out 2>/dev/null | head -1)
  [ -z "$LOG" ] && LOG=$(ls -t $CODE/logs/*_${JID}.out 2>/dev/null | head -1)
  if [ -z "$LOG" ] || [ ! -f "$LOG" ]; then echo "$JID $JN $ELAP $NODE log=NONE"; continue; fi
  AGE=$(( ( $(date +%s) - $(stat -c %Y "$LOG") ) / 60 ))
  CH=$(tr '\r' '\n' < "$LOG" | grep -o 'Chunk-[0-9]*/Total-chunk-[0-9]*' | tail -1)
  echo "$JID $JN elapsed=$ELAP node=$NODE log_age_min=$AGE ${CH:-chunk=none}"
done

echo "--- FAILED / non-zero exits in last 400 finished jobs ---"
sacct -u "$USER" -S $(date -d '2 days ago' +%F) -X \
  --format=JobID%12,JobName%34,State%14,Elapsed,MaxRSS,NodeList%14 2>/dev/null \
  | grep -E 'condunit|u_functional' | grep -vE 'COMPLETED|RUNNING|PENDING' || echo "(none)"

echo "--- CONDITIONING AUDIT (per completed run log) ---"
NCONF=0; NBAD=0
for S in $(cat $CODE/samples.txt); do
  for L in "$WORKDIR/$S/cond_${ARM}"_*_run.log; do
    [ -e "$L" ] || continue
    if grep -Fq "Baseline annotation is not provided" "$L"; then echo "UNCONDITIONED: $L"; NBAD=$((NBAD+1))
    elif grep -Fq "Using additional baseline annotations" "$L"; then NCONF=$((NCONF+1))
    else echo "NO-CONFIRM (may be mid-run): $L"; fi
  done
done
echo "conditioning confirmed: $NCONF   unconditioned: $NBAD"

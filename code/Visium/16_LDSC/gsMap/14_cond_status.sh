#!/bin/bash
# 14_cond_status.sh <ARM> -- terse status: counts, anomalies only.
# Prints a running unit ONLY if it looks unhealthy (log older than 45 min, or no
# chunk progress after 20 min of elapsed time). Healthy units are summarised.
set -uo pipefail
ARM="${1:-functional}"
CODE=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
WORKDIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}

NOUT=$(while read -r S; do ls -1 "$WORKDIR/$S/spatial_ldsc/"*.csv.gz 2>/dev/null; done < $CODE/samples.txt | wc -l)
NPEND=0
while read -r T; do while read -r S; do
  [ -s "$WORKDIR/$S/spatial_ldsc/${S}_${T}.csv.gz" ] || NPEND=$((NPEND+1))
done < $CODE/samples.txt; done < $CODE/traits.txt
echo "OUTPUTS: $NOUT/280   still-missing: $NPEND"
echo "unit jobs: running=$(squeue -u $USER -h -t R -o %j | grep -c "^u_${ARM}_") pending=$(squeue -u $USER -h -t PD -o %j | grep -c "^u_${ARM}_")"

NSLOW=0
squeue -u "$USER" -h -t R -o "%i %j %M %N" | grep " u_${ARM}_" | while read -r JID JN ELAP NODE; do
  LOG=$(ls -t $CODE/logs/${JN}_${JID}.out 2>/dev/null | head -1)
  [ -f "$LOG" ] || { echo "ANOMALY $JID $JN node=$NODE elapsed=$ELAP log=NONE"; continue; }
  AGE=$(( ( $(date +%s) - $(stat -c %Y "$LOG") ) / 60 ))
  CH=$(tr '\r' '\n' < "$LOG" | grep -o 'Chunk-[0-9]*/Total-chunk-[0-9]*' | tail -1)
  if [ "$AGE" -gt 45 ]; then echo "STALL $JID $JN node=$NODE elapsed=$ELAP log_age_min=$AGE ${CH:-chunk=none}"; fi
done

echo "--- non-COMPLETED unit jobs (last 2 days) ---"
sacct -u "$USER" -S $(date -d '2 days ago' +%F) -X --format=JobID%12,JobName%32,State%16,Elapsed,MaxRSS,NodeList%14 2>/dev/null \
  | grep "u_${ARM}_" | grep -vE 'COMPLETED|RUNNING|PENDING' || echo "(none)"

echo "--- conditioning audit ---"
NCONF=0; NBAD=0; NNC=0
for S in $(cat $CODE/samples.txt); do
  for L in "$WORKDIR/$S/cond_${ARM}"_*_run.log; do
    [ -e "$L" ] || continue
    if grep -Fq "Baseline annotation is not provided" "$L"; then echo "UNCONDITIONED: $L"; NBAD=$((NBAD+1))
    elif grep -Fq "Using additional baseline annotations" "$L"; then NCONF=$((NCONF+1))
    else NNC=$((NNC+1)); fi
  done
done
echo "confirmed=$NCONF unconditioned=$NBAD no-confirm-yet(mid-run)=$NNC"

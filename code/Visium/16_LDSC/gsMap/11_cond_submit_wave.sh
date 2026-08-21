#!/bin/bash
# 11_cond_submit_wave.sh -- top the queue up to a target depth with pending units.
#
#   bash 11_cond_submit_wave.sh <ARM> <QUEUE_TARGET> [EXTRA_EXCLUDE]
#
# Fairshare on this account is 0.031 with 300+ jobs typically pending and a
# sibling track running the other arm, so we hold our own queue at <= QUEUE_TARGET
# rather than dumping all 280 units in at once.
#
# A unit is PENDING if its per-spot output is absent/empty AND no job for it is
# already queued or running. Job names encode the unit, so squeue acts as the
# lock; cond_skip_units.txt covers units in flight under a legacy job name that
# the name-lock cannot see. Safe to re-run: it never double-submits.
set -uo pipefail
ARM="${1:?usage: 11_cond_submit_wave.sh ARM QUEUE_TARGET [EXTRA_EXCLUDE]}"
TARGET="${2:?}"
EXTRA="${3:-}"

CODE=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
WORKDIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}
EXCL="compute-058,compute-175,compute-053"
[ -n "$EXTRA" ] && EXCL="$EXCL,$EXTRA"

mkdir -p "$CODE/logs"
INFLIGHT=$(squeue -u "$USER" -h -o "%j" | tr '\n' ' ')

# pilot traits first (they have unconditioned counterparts to compare against),
# then traits.txt order; samples in samples.txt order within each trait.
ORDER=$( { echo MDD; echo SCZ; echo Height; grep -vxE 'MDD|SCZ|Height' "$CODE/traits.txt"; } )

# Count only THIS arm's unit jobs toward the target. A sibling track runs the
# other arm under the same account; counting its jobs here would starve this
# sweep whenever the sibling is busy. The cap is "<= ~60 of my own".
NQ=$(squeue -u "$USER" -h -o "%j" | grep -c "^u_${ARM}_")
echo "queue depth now: $NQ   target: $TARGET   exclude: $EXCL"
nsub=0
for TRAIT in $ORDER; do
  while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue
    [ "$NQ" -ge "$TARGET" ] && { echo "queue target reached; submitted $nsub this wave"; exit 0; }
    SPOT="$WORKDIR/$SAMPLE/spatial_ldsc/${SAMPLE}_${TRAIT}.csv.gz"
    [ -s "$SPOT" ] && continue
    JN="u_${ARM}_${SAMPLE}_${TRAIT}"
    case " $INFLIGHT " in *" $JN "*) continue ;; esac
    if [ -f "$CODE/cond_skip_units.txt" ] && grep -qxF "$SAMPLE $TRAIT" "$CODE/cond_skip_units.txt"; then
      echo "skip (external job in flight): $SAMPLE $TRAIT"; continue; fi
    JID=$(sbatch --parsable \
      --job-name="$JN" \
      --exclude="$EXCL" \
      --output="$CODE/logs/${JN}_%j.out" \
      --error="$CODE/logs/${JN}_%j.out" \
      --export=ALL,ARM="$ARM",SAMPLE="$SAMPLE",TRAIT="$TRAIT" \
      "$CODE/10_cond_unit.sh")
    echo "submitted $JID $JN"
    nsub=$((nsub+1)); NQ=$((NQ+1))
  done < "$CODE/samples.txt"
done
echo "no pending units left to submit; submitted $nsub this wave"

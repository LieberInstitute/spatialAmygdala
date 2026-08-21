#!/bin/bash
#SBATCH --job-name=scdrs_arr
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --array=1-42%20
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/06b_score.%A_%a.log
# -----------------------------------------------------------------------------
# 06b_score_array.sh -- FULL SWEEP as a per-trait array. Alternative to 06.
#
#   sbatch 06b_score_array.sh              # all 42 traits, 20 at a time
#   sbatch --array=35 06b_score_array.sh   # just one trait (task = .gs row)
#
# WHY AN ARRAY RATHER THAN THE SINGLE JOB IN 06
#   06's header comment argued that one job amortises the matrix load across
#   traits. Measured on job 34721738 that argument does not hold:
#
#     setup (read 224k x 19k h5ad, normalise, match genes) : ~7 min, ONCE
#     per-trait control-set sampling (1000 sets)           : ~32-39 min
#
#   Setup is ~15-20% of a single trait's cost, so paying it 42 times adds
#   ~5 h of CPU but removes the serial chain. 42 x ~35 min serial is ~25 h;
#   the same work 20-wide is ~3 waves, ~2 h wall-clock.
#
#   Measured on the same job, scdrs compute-score is effectively
#   SINGLE-THREADED (AveCPU 2:08:48 against 2:10:15 elapsed = ~1.0 core with
#   8 allocated), so cpus-per-task drops 8 -> 2 and 20 concurrent tasks cost
#   ~20 cores of the 3592-CPU partition rather than 160.
#
#   MEMORY. Peak RSS scales weakly with trait count because the 224k x 19k
#   sparse matrix dominates and the gene sets are small:
#
#     3 traits  (pilot 34715015) : 14.1 GB
#     42 traits (job   34721738) : 22.0 GB
#     => ~13.5 GB baseline + ~200 MB per trait
#
#   A ONE-TRAIT task therefore peaks near 13.7 GB, so 24 G gives ~1.75x
#   headroom. Note this is per task: 20 concurrent tasks reserve ~480 G across
#   the partition, versus 32 G for the single job -- the wall-clock win is paid
#   for in reserved memory, which is why 24 G and not 32 G.
#
# IDEMPOTENT: a task whose $trait.score.gz already exists exits immediately,
# so this can be submitted after a partial run of 06 and will only fill gaps.
# Nothing is ever deleted.
# -----------------------------------------------------------------------------
set -euo pipefail

# $0 does NOT resolve to this file under sbatch (SLURM runs a spool copy), so
# the code directory is hardcoded. Override with CODE_DIR=... if relocated.
CODE_DIR=${CODE_DIR:-/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS}
source "$CODE_DIR/config.sh"

echo "**** Job starts ****"; date; hostname
activate_env
mkdir -p "$SCORE_DIR" "$CODE/logs" "$SCRATCH/gs_split"

GS=${GS_CHUNK:-$GS_DIR/magma_top${GS_NMAX}.gs}

# .gs is  header + one row per trait, so data row N is file line N+1.
LINE=$(( ${SLURM_ARRAY_TASK_ID:-1} + 1 ))
TRAIT=$(sed -n "${LINE}p" "$GS" | cut -f1)

if [ -z "${TRAIT:-}" ]; then
    echo "no trait at .gs line $LINE -- array index out of range"; exit 0
fi
echo "task ${SLURM_ARRAY_TASK_ID:-1} -> trait '$TRAIT' (.gs line $LINE)"

# Skip only if BOTH outputs exist AND both are complete gzip streams. Testing
# `-s` alone is not enough: a job killed mid-write leaves a large but truncated
# .gz that a size check accepts and the next step then fails on with
# "EOFError: Compressed file ended before the end-of-stream marker was
# reached". That happened to AUD_EA_MVP when the serial sweep was cancelled.
SKIP=1
for F in "$SCORE_DIR/$TRAIT.score.gz" "$SCORE_DIR/$TRAIT.full_score.gz"; do
    if [ ! -s "$F" ] || ! gzip -t "$F" 2>/dev/null; then SKIP=0; fi
done
if [ "$SKIP" = "1" ]; then
    echo "$TRAIT already scored and both files pass gzip -t -- nothing to do"
    ls -lh "$SCORE_DIR/$TRAIT".*score.gz
    echo "**** Job ends ****"; date; exit 0
fi
echo "$TRAIT needs (re)scoring -- missing or truncated output"

# Single-trait .gs for this task. Written to scratch and left in place.
ONE=$SCRATCH/gs_split/${TRAIT}.gs
head -1 "$GS" > "$ONE"
sed -n "${LINE}p" "$GS" >> "$ONE"
echo "wrote $ONE ($(wc -l < "$ONE") lines)"

/usr/bin/time -v "$SCDRS" compute-score \
    --h5ad-file "$H5AD_SCDRS" \
    --h5ad-species human \
    --gs-file "$ONE" \
    --gs-species human \
    --cov-file "$COV_FILE" \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl "$N_CTRL" \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    --out-folder "$SCORE_DIR"

ls -lh "$SCORE_DIR/$TRAIT".*
echo "**** Job ends ****"; date

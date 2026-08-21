#!/bin/bash
#SBATCH --job-name=cond_neuro_cauchy
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_neuro_cauchy_%j.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_neuro_cauchy_%j.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --exclude=compute-058,compute-175,compute-053
# -----------------------------------------------------------------------------
# 11_cond_neuronal_cauchy_across.sh -- pool spot-level p-values across all 7
# capture areas into ONE p-value per (domain, trait) for ARM=neuronal (test 3).
# Mirrors 05_cauchy_across_samples.sh but reads the conditional workdir.
#
# Idempotent: a trait is skipped when its output already exists, and a trait
# whose 7 per-sample spatial_ldsc outputs are not all present is DEFERRED, not
# failed -- so this can be re-run after each sweep wave.
# -----------------------------------------------------------------------------
set -uo pipefail

PROJ=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala
CODE=$PROJ/code/Visium/16_LDSC/gsMap
ARM=neuronal
WORKDIR=$PROJ/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}
SCRATCH=/users/mtotty/claude_scratch/gsmap
ANNOT=BS_k16_Semisupervised_wAI
OUTDIR=$WORKDIR/cauchy_across_samples
mkdir -p "$OUTDIR" "$CODE/logs"

echo "**** Job starts ****"; date
echo "host: $(hostname)"
module load conda/3-24.3.0
source activate $SCRATCH/envs/gsmap
echo "gsmap: $(which gsmap)"

SAMPLE_LIST=$(tr '\n' ' ' < $CODE/samples.txt)
echo "samples: $SAMPLE_LIST"

n_ok=0; n_fail=0; n_skip=0; n_defer=0
while read -r TRAIT; do
    [ -z "$TRAIT" ] && continue
    OUT=$OUTDIR/${TRAIT}_cauchy.csv.gz
    if [ -s "$OUT" ]; then echo "SKIP $TRAIT (present)"; n_skip=$((n_skip+1)); continue; fi
    nhave=0
    for S in $SAMPLE_LIST; do
        [ -s "$WORKDIR/$S/spatial_ldsc/${S}_${TRAIT}.csv.gz" ] && nhave=$((nhave+1))
    done
    if [ "$nhave" -ne 7 ]; then
        echo "DEFER $TRAIT ($nhave/7 samples ready)"; n_defer=$((n_defer+1)); continue
    fi
    echo "==== $TRAIT ===="; date
    gsmap run_cauchy_combination \
        --workdir "$WORKDIR" \
        --trait_name "$TRAIT" \
        --annotation "$ANNOT" \
        --sample_name_list $SAMPLE_LIST \
        --output_file "$OUT"
    RC=$?
    if [ -s "$OUT" ]; then echo "OK   $TRAIT -> $OUT (rc=$RC)"; n_ok=$((n_ok+1))
    else echo "FAIL $TRAIT (rc=$RC)"; n_fail=$((n_fail+1)); fi
done < $CODE/traits.txt

echo "----"
echo "ok=$n_ok fail=$n_fail skipped=$n_skip deferred=$n_defer"
ls -la "$OUTDIR" | tail -50
echo "**** Job ends ****"; date

#!/bin/bash
#SBATCH --job-name=gsmap_cond_trait
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_trait_%A_%a.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_trait_%A_%a.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-3
set -euo pipefail

ARM="${ARM:?set ARM=functional or ARM=amygdala}"
SAMPLE="${SAMPLE:-Br6471}"
WORKDIR="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}"
export PATH="/users/mtotty/claude_scratch/gsmap/gsMap_resource/../envs/gsmap/bin:$PATH"

# pilot traits: ones we already have UNCONDITIONED results for
TRAIT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/cond_pilot_traits.txt)
GWAS="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap/GWAS/$(sed -n "${SLURM_ARRAY_TASK_ID}p" /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/cond_pilot_files.txt)"

echo "ARM=$ARM TRAIT=$TRAIT SAMPLE=$SAMPLE"
echo "--- verifying the conditioning is ACTUALLY active ---"
ls -d "$WORKDIR/$SAMPLE/generate_ldscore/additional_baseline" || {
  echo "FATAL: additional_baseline dir absent -> gsMap would SILENTLY run unconditioned"; exit 4; }

RUNLOG="$WORKDIR/$SAMPLE/cond_${ARM}_${TRAIT}_run.log"
mkdir -p "$(dirname "$RUNLOG")"

# Run to completion, capturing the full log. No pipe into grep: a short-circuiting
# reader (grep -q) would SIGPIPE gsmap mid-run and, under pipefail, abort the job.
set +e   # capture rc ourselves; `|| true` would reset PIPESTATUS to 0
gsmap run_spatial_ldsc \
  --workdir "$WORKDIR" \
  --sample_name "$SAMPLE" \
  --trait_name "$TRAIT" \
  --sumstats_file "$GWAS" \
  --num_processes 8 \
  --use_additional_baseline_annotation True 2>&1 | tee "$RUNLOG"
LDSC_RC=${PIPESTATUS[0]}
set -e
echo "gsmap run_spatial_ldsc exited rc=$LDSC_RC (non-fatal; success is judged by outputs below)"

# Verify AFTER the run finishes, against the saved log.
if grep -Fq "Using additional baseline annotations" "$RUNLOG"; then
  echo "CONFIRMED: conditioning active"
else
  echo "FATAL: log never said 'Using additional baseline annotations' -> ran UNCONDITIONED"
  exit 5
fi

# gsMap can exit non-zero on cosmetic post-processing (the pilot died in report
# generation with all science complete). Never let rc alone decide.
SPOT="$WORKDIR/$SAMPLE/spatial_ldsc/${SAMPLE}_${TRAIT}.csv.gz"
if [ ! -s "$SPOT" ]; then
  echo "FAIL: no per-spot output $SPOT (rc=$LDSC_RC)"; exit 6
fi
echo "per-spot output present: $SPOT"

CAU_RC=0
gsmap run_cauchy_combination \
  --workdir "$WORKDIR" --sample_name "$SAMPLE" \
  --trait_name "$TRAIT" --annotation BS_k16_Semisupervised_wAI || CAU_RC=$?

OUT="$WORKDIR/$SAMPLE/cauchy_combination/${SAMPLE}_${TRAIT}.Cauchy.csv.gz"
if [ -s "$OUT" ]; then
  echo "OK $TRAIT -> $OUT (ldsc rc=$LDSC_RC cauchy rc=$CAU_RC)"
else
  echo "FAIL: no Cauchy output for $TRAIT (cauchy rc=$CAU_RC)"; exit 6
fi

#!/bin/bash
#SBATCH --job-name=cond_multi
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_multi_%j.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_multi_%j.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --exclude=compute-053,compute-058,compute-175,compute-142,compute-148,compute-103
# -----------------------------------------------------------------------------
# PROTOTYPE: run MANY traits for ONE sample in a single gsMap invocation.
# gsMap's chunk loop iterates traits INSIDE each chunk, so the per-chunk
# LD-score matmul is computed once and reused across all traits in the config.
# Requires: ARM, SAMPLE, CFG (yaml of trait -> sumstats path)
# -----------------------------------------------------------------------------
set -uo pipefail
PROJ=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala
CODE=$PROJ/code/Visium/16_LDSC/gsMap
ARM="${ARM:?set ARM=functional or ARM=neuronal}"
SAMPLE="${SAMPLE:?set SAMPLE}"
CFG="${CFG:?set CFG=path to sumstats config yaml}"
WORKDIR=$PROJ/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}
export PATH="/users/mtotty/claude_scratch/gsmap/envs/gsmap/bin:$PATH"

case "$ARM" in
  functional) ANNOT=/users/mtotty/claude_scratch/gsmap/additional_annotation/gsMap_additional_annotation ;;
  neuronal)   ANNOT=$PROJ/processed-data/Visium/16_LDSC/gsMap_cond/annotation_neuronal ;;
  *) echo "FATAL: bad ARM=$ARM"; exit 3 ;;
esac

LD=$WORKDIR/$SAMPLE/generate_ldscore/additional_baseline
[ -d "$LD" ] || { echo "FATAL: missing $LD -- would run UNCONDITIONED"; exit 4; }
NF=$(ls "$LD"/baseline.*.l2.ldscore.feather 2>/dev/null | wc -l)
[ "$NF" -eq 22 ] || { echo "FATAL: $NF/22 baseline feathers"; exit 5; }
[ -e "$WORKDIR/$SAMPLE/find_latent_representations/${SAMPLE}_add_latent.h5ad" ] \
  || { echo "FATAL: add_latent.h5ad missing (cauchy would fail)"; exit 7; }

NTR=$(grep -cE '^[A-Za-z0-9_]+:' "$CFG")
echo "=== multi-trait run: ARM=$ARM SAMPLE=$SAMPLE traits=$NTR ==="
date '+start %F %T'
T0=$(date +%s)

RUNLOG=$CODE/logs/cond_multi_${SAMPLE}_${ARM}_${SLURM_JOB_ID}.gsmap.log
set +e
gsmap run_spatial_ldsc \
  --workdir "$WORKDIR" --sample_name "$SAMPLE" \
  --sumstats_config_file "$CFG" \
  --additional_baseline_annotation "$ANNOT" \
  --num_processes 8 2>&1 | tee "$RUNLOG"
LDSC_RC=${PIPESTATUS[0]}
set -e
echo "gsmap run_spatial_ldsc rc=$LDSC_RC (success judged by outputs)"

grep -q "Using additional baseline annotation" "$RUNLOG" \
  || echo "WARNING: no 'Using additional baseline annotation' line -- CHECK CONDITIONING"
grep -q "Baseline annotation is not provided" "$RUNLOG" \
  && { echo "FATAL: ran UNCONDITIONED"; exit 6; }

T1=$(date +%s)
echo "=== spatial_ldsc wall: $(( (T1-T0)/60 )) min for $NTR traits ==="

NOUT=0
while IFS= read -r TR; do
  OUT=$WORKDIR/$SAMPLE/spatial_ldsc/${SAMPLE}_${TR}.csv.gz
  if [ -s "$OUT" ]; then
    NOUT=$((NOUT+1))
    gsmap run_cauchy_combination --workdir "$WORKDIR" --sample_name "$SAMPLE" \
      --trait_name "$TR" --annotation BS_k16_Semisupervised_wAI >/dev/null 2>&1 \
      || echo "  cauchy failed: $TR"
  else
    echo "  MISSING output: $TR"
  fi
done < <(grep -oE '^[A-Za-z0-9_]+:' "$CFG" | tr -d ':')

T2=$(date +%s)
echo "=== outputs: $NOUT/$NTR   total wall: $(( (T2-T0)/60 )) min ==="
echo "=== per-trait cost: $(( (T2-T0)/60/NTR )) min/trait (vs ~127 min single-trait) ==="
date '+end %F %T'

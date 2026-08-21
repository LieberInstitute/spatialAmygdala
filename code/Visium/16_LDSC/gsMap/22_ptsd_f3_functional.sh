#!/bin/bash
#SBATCH --job-name=ptsd_f3
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/ptsdf3_%A_%a.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/ptsdf3_%A_%a.err
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=36G
#SBATCH --time=12:00:00
#SBATCH --array=1-7
#SBATCH --exclude=compute-053,compute-058,compute-175
#
# 22_ptsd_f3_functional.sh -- PTSD_F3 under the FUNCTIONAL-conditioned arm (test 2).
# One task per donor. Same pipeline as 10_cond_unit.sh; only the trait differs.
#
# RESOURCES -- 36G/12h, NOT the 48G/1-day the arm was originally submitted at. That
# right-sizing was applied to 10_cond_unit.sh and 12_cond_neuronal_unit.sh mid-sweep
# (backups kept as *.bak48G) because 48G was scheduling badly: the `shared` partition
# bills CPU=1.0 and Mem=1.0/MB, so billing was 99.98% memory and 105 units sat pending
# behind it. Dropping to 36G cut billing 49,160 -> 36,872 (-25%) and several units
# scheduled immediately. Basis: observed MaxRSS 23-27 GB (~33% headroom at 36G; 24G
# OOM-killed every pilot task), median unit 2.11 h and p90 3.90 h against a 12 h cap,
# median CPU peak 634% across 91 jobs so 8 cores stay warranted.
# If a task OOMs or hits the wall, raise to 48G/1-00:00:00 rather than retrying as-is.
# compute-053 runs ~half rate; 058/175 stall silently.
# SHELL FLAGS: no `set -e`, no `pipefail`. gsMap can exit non-zero having written correct
# output (cosmetic step-6 report bug) and can exit ZERO having written nothing, so an
# aborting shell would skip the output checks that actually decide success. Every guard
# below therefore tests explicitly and calls `exit 1` itself.

G="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap"
PD="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC"
WORKDIR="$PD/gsMap_cond_functional"
TRAIT="PTSD_F3"
GWAS="$PD/gsMap/GWAS/${TRAIT}.gz"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$G/samples.txt")

export PATH="/users/mtotty/claude_scratch/gsmap/envs/gsmap/bin:$PATH"
echo "sample=$SAMPLE trait=$TRAIT host=$(hostname) start=$(date)"

[ -f "$GWAS" ] || { echo "FATAL: $GWAS missing -- run 21_convert_ptsd_f3.sh first" >&2; exit 1; }

# Both prerequisites gsMap reads but does NOT create in a conditional workdir.
# Missing either kills the run -- the Cauchy one only AFTER the expensive regression.
for f in "latent_to_gene/${SAMPLE}_gene_marker_score.feather" \
         "find_latent_representations/${SAMPLE}_add_latent.h5ad"; do
  [ -e "$WORKDIR/$SAMPLE/$f" ] || { echo "FATAL: missing prerequisite $f" >&2; exit 1; }
done

# Guard the silent-disable trap: config.py turns conditioning OFF without error when
# additional_baseline/ is absent, exits 0, and writes UNCONDITIONED results.
AB="$WORKDIR/$SAMPLE/generate_ldscore/additional_baseline"
n_ab=$(ls "$AB"/baseline.*.l2.ldscore.feather 2>/dev/null | wc -l)
[ "$n_ab" -eq 22 ] || { echo "FATAL: expected 22 additional_baseline feathers, found $n_ab" >&2; exit 1; }

RUNLOG="$G/logs/ptsdf3_${SAMPLE}_${SLURM_JOB_ID}.gsmap.log"
gsmap run_spatial_ldsc \
  --workdir "$WORKDIR" \
  --sample_name "$SAMPLE" \
  --trait_name "$TRAIT" \
  --sumstats_file "$GWAS" \
  --num_processes 8 \
  --use_additional_baseline_annotation True 2>&1 | tee "$RUNLOG"
LDSC_RC=${PIPESTATUS[0]}   # must be read on the line immediately after the pipeline
echo "run_spatial_ldsc rc=$LDSC_RC (non-fatal; success judged by outputs)"

if grep -q "Baseline annotation is not provided" "$RUNLOG"; then
  echo "FATAL: ran UNCONDITIONED -- conditioning silently disabled" >&2; exit 1
fi

OUT="$WORKDIR/$SAMPLE/spatial_ldsc/${SAMPLE}_${TRAIT}.csv.gz"
[ -s "$OUT" ] || { echo "FATAL: no per-spot output at $OUT" >&2; exit 1; }

gsmap run_cauchy_combination \
  --workdir "$WORKDIR" --sample_name "$SAMPLE" \
  --trait_name "$TRAIT" --annotation BS_k16_Semisupervised_wAI

echo "done=$(date) out=$(du -h "$OUT" | cut -f1)"

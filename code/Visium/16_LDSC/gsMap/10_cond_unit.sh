#!/bin/bash
#SBATCH --job-name=condunit
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=36G
#SBATCH --time=12:00:00
#SBATCH --exclude=compute-058,compute-175,compute-053
# -----------------------------------------------------------------------------
# 10_cond_unit.sh -- ONE (trait, sample) conditional spatial-LDSC unit.
#
# Derived from 09_cond_spatial_ldsc.sh, which took TRAIT from a SLURM array index
# into cond_pilot_traits.txt. Here TRAIT arrives via --export so any of the 40
# traits in traits.txt can be run, one unit per job.
#
#   sbatch --export=ALL,ARM=functional,SAMPLE=Br6660,TRAIT=SCZ \
#          --output=$CODE/logs/unit_functional_Br6660_SCZ_%j.out \
#          $CODE/10_cond_unit.sh
#
# mem: siblings measured 27.5G / 24.5G MaxRSS via sacct at 8 cpus; 24G was killed
# OUT_OF_MEMORY. gsMap's own "Peak memory usage" line under-reports MaxRSS by
# 1.55-1.78x, so 48G is sized off sacct, not off the gsMap log.
#
# Exit codes: 3 = GWAS path not resolvable, 4 = additional_baseline dir missing
# (gsMap would SILENTLY run unconditioned), 5 = log never confirmed conditioning,
# 6 = per-spot output absent. Success is judged by OUTPUT FILE, never by gsMap's
# own rc: it exits 1 on cosmetic report errors with all science complete, and
# exits 0 having computed nothing.
# -----------------------------------------------------------------------------
set -uo pipefail

ARM="${ARM:?set ARM=functional or ARM=neuronal}"
SAMPLE="${SAMPLE:?set SAMPLE}"
TRAIT="${TRAIT:?set TRAIT}"

CODE=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
WORKDIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}
export PATH="/users/mtotty/claude_scratch/gsmap/envs/gsmap/bin:$PATH"

echo "**** unit start ****"; date
echo "ARM=$ARM SAMPLE=$SAMPLE TRAIT=$TRAIT host=$(hostname) job=${SLURM_JOB_ID:-NA}"
echo "gsmap: $(which gsmap)"

SPOT="$WORKDIR/$SAMPLE/spatial_ldsc/${SAMPLE}_${TRAIT}.csv.gz"

# Idempotent: a requeued unit whose output already landed does no work.
if [ -s "$SPOT" ]; then
  echo "SKIP: per-spot output already present and non-empty: $SPOT"
  ls -l "$SPOT"
  exit 0
fi

# GWAS path from the single source of truth, gwas_config_full.yaml.
GWAS=$(awk -v t="$TRAIT" '$1==t":" {print $2}' "$CODE/gwas_config_full.yaml")
if [ -z "$GWAS" ] || [ ! -s "$GWAS" ]; then
  echo "FATAL: no readable GWAS for TRAIT=$TRAIT (resolved: '$GWAS')"; exit 3; fi
echo "GWAS=$GWAS"

echo "--- verifying the conditioning is ACTUALLY active ---"
ls -d "$WORKDIR/$SAMPLE/generate_ldscore/additional_baseline" || {
  echo "FATAL: additional_baseline dir absent -> gsMap would SILENTLY run unconditioned"; exit 4; }

# Pre-flight for the CAUCHY step, checked BEFORE the 2.5-3h regression:
# run_cauchy_combination reads <S>_add_latent.h5ad from THIS workdir, which the
# conditional workdirs did not originally carry. It is symlinked from arm 0.
# Failing fast here costs seconds; discovering it after the regression wastes hours.
LATENT="$WORKDIR/$SAMPLE/find_latent_representations/${SAMPLE}_add_latent.h5ad"
if [ ! -s "$LATENT" ] || [ ! -r "$LATENT" ]; then
  echo "FATAL: latent h5ad missing/unreadable -> cauchy step would die: $LATENT"; exit 9; fi
echo "latent h5ad present: $(readlink -f "$LATENT")"

RUNLOG="$WORKDIR/$SAMPLE/cond_${ARM}_${TRAIT}_run.log"
mkdir -p "$(dirname "$RUNLOG")"

# No pipe into grep: a short-circuiting reader would SIGPIPE gsmap mid-run.
gsmap run_spatial_ldsc \
  --workdir "$WORKDIR" \
  --sample_name "$SAMPLE" \
  --trait_name "$TRAIT" \
  --sumstats_file "$GWAS" \
  --num_processes 8 \
  --use_additional_baseline_annotation True 2>&1 | tee "$RUNLOG"
LDSC_RC=${PIPESTATUS[0]}
echo "gsmap run_spatial_ldsc exited rc=$LDSC_RC (non-fatal; success judged by outputs)"

# Verify AFTER the run, against the saved log -- both directions.
if grep -Fq "Baseline annotation is not provided" "$RUNLOG"; then
  echo "FATAL: log says 'Baseline annotation is not provided' -> ran UNCONDITIONED"; exit 5; fi
if grep -Fq "Using additional baseline annotations" "$RUNLOG"; then
  echo "CONFIRMED: conditioning active"
else
  echo "FATAL: log never said 'Using additional baseline annotations' -> ran UNCONDITIONED"; exit 5; fi

if [ ! -s "$SPOT" ]; then
  echo "FAIL: no per-spot output $SPOT (rc=$LDSC_RC)"; exit 6; fi
echo "per-spot output present:"; ls -l "$SPOT"

# Per-sample Cauchy is a convenience; the deliverable is the cross-sample pool.
# A cosmetic failure here must NOT invalidate a good per-spot result.
CAU_RC=0
gsmap run_cauchy_combination \
  --workdir "$WORKDIR" --sample_name "$SAMPLE" \
  --trait_name "$TRAIT" --annotation BS_k16_Semisupervised_wAI || CAU_RC=$?
CAU="$WORKDIR/$SAMPLE/cauchy_combination/${SAMPLE}_${TRAIT}.Cauchy.csv.gz"
if [ -s "$CAU" ]; then echo "per-sample cauchy OK -> $CAU"
else echo "WARN: per-sample cauchy missing (rc=$CAU_RC) -- per-spot result stands"; fi

echo "UNIT_OK $ARM $SAMPLE $TRAIT"
echo "**** unit end ****"; date

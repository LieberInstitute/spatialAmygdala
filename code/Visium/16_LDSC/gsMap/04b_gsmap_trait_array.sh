#!/bin/bash
#SBATCH --job-name=gsmap_trait
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/trait_%A_%a.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/trait_%A_%a.err
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-40
# -----------------------------------------------------------------------------
# 04b_gsmap_trait_array.sh -- STAGE 2 of 2: the full sweep.
# One array task per TRAIT (40 tasks); each runs all 7 capture areas.
#
# WHY quick_mode AND NOT run_spatial_ldsc
# ---------------------------------------
# The first version of this script called `gsmap run_spatial_ldsc` directly.
# That is WRONG for a quick_mode workdir and fails SILENTLY:
#   run_spatial_ldsc expects CHUNKED ldscore files in {sample}/generate_ldscore/.
#   In quick_mode that directory holds only symlinks to the pre-built resource
#   (baseline/, SNP_gene_pair/) plus a .done marker -- no chunks. So it logs
#   "Find 0 chunked files ... Spatial LDSC finished!", EXITS 0, and writes nothing.
#   Job 34696332 completed 37/40 tasks in ~90 s each having computed nothing.
#
# `gsmap quick_mode` is cache-aware: with steps 1-2 already built by
# 04a_gsmap_cache_pass.sh it reports "Step 1/2 completed in 0h 0m" and does only
# the step-4 regression. Verified in the Br6471 07:20 pipeline log.
#
# EXIT CODES ARE NOT TRUSTWORTHY HERE
# -----------------------------------
# quick_mode's final report step (step 6) crashes with KeyError: 'Z' in
# gsMap/diagnosis.py -- it re-reads the raw sumstats expecting an uppercase Z
# column. This is COSMETIC: spatial_ldsc (step 4) and cauchy (step 5) have
# already written their outputs. So success is judged by the OUTPUT FILE
# existing, never by the exit code.
#
#   sbatch 04b_gsmap_trait_array.sh
#   sbatch --array=34 04b_gsmap_trait_array.sh     # re-run just SCZ
# -----------------------------------------------------------------------------
set -uo pipefail

PROJ=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala
CODE=$PROJ/code/Visium/16_LDSC/gsMap
WORKDIR=$PROJ/processed-data/Visium/16_LDSC/gsMap
SCRATCH=/users/mtotty/claude_scratch/gsmap

ANNOT=BS_k16_Semisupervised_wAI
CONFIG=$CODE/gwas_config_full.yaml
SLICE_MEAN=$WORKDIR/sample_slice_mean.parquet

TRAIT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CODE/traits.txt)
if [ -z "$TRAIT" ]; then
    echo "ERROR: no trait at line ${SLURM_ARRAY_TASK_ID} of traits.txt"; exit 1
fi
SUMSTATS=$(grep -E "^${TRAIT}:" $CONFIG | head -1 | sed 's/^[^:]*:[[:space:]]*//')
if [ -z "$SUMSTATS" ] || [ ! -e "$SUMSTATS" ]; then
    echo "ERROR: sumstats for $TRAIT not found (resolved to '$SUMSTATS')"; exit 1
fi
if [ ! -f "$SLICE_MEAN" ]; then
    echo "ERROR: $SLICE_MEAN missing -- run 03_create_slice_mean.sh first"; exit 1
fi

mkdir -p "$CODE/logs" "$WORKDIR/cauchy_across_samples"
echo "**** Job starts ****"; date
echo "trait: $TRAIT"; echo "sumstats: $SUMSTATS"
echo "host: $(hostname)  cpus: ${SLURM_CPUS_PER_TASK:-NA}"

module load conda/3-24.3.0
source activate $SCRATCH/envs/gsmap

MISSING=0
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue
    OUT=$WORKDIR/$SAMPLE/spatial_ldsc/${SAMPLE}_${TRAIT}.csv.gz
    if [ -s "$OUT" ]; then
        echo "--- $SAMPLE / $TRAIT already present, skipping"; continue
    fi
    if [ ! -d "$WORKDIR/$SAMPLE/latent_to_gene" ]; then
        echo "ERROR: $SAMPLE has no cached latent_to_gene -- run 04a first"
        MISSING=1; continue
    fi
    echo "--- quick_mode (cached steps 1-2): $SAMPLE / $TRAIT"; date
    gsmap quick_mode \
        --workdir "$WORKDIR" \
        --sample_name "$SAMPLE" \
        --gsMap_resource_dir "$SCRATCH/gsMap_resource" \
        --hdf5_path "$WORKDIR/ST/${SAMPLE}.h5ad" \
        --annotation "$ANNOT" \
        --data_layer count \
        --sumstats_file "$SUMSTATS" \
        --trait_name "$TRAIT" \
        --gM_slices "$SLICE_MEAN" \
        --max_processes ${SLURM_CPUS_PER_TASK:-8}
    # exit code ignored on purpose: step 6 (report) crashes cosmetically.
    if [ -s "$OUT" ]; then
        echo "    OK  $(stat -c %s "$OUT") bytes"
    else
        echo "    FAILED: $OUT not written"; MISSING=1
    fi
done < $CODE/samples.txt

if [ "$MISSING" -eq 0 ]; then
    echo "--- cauchy combination across samples: $TRAIT"; date
    gsmap run_cauchy_combination \
        --workdir "$WORKDIR" \
        --trait_name "$TRAIT" \
        --annotation "$ANNOT" \
        --sample_name_list $(tr '\n' ' ' < $CODE/samples.txt) \
        --output_file "$WORKDIR/cauchy_across_samples/${TRAIT}_cauchy.csv.gz" \
        || echo "ERROR: cauchy combination failed for $TRAIT"
else
    echo "SKIPPING cauchy for $TRAIT -- one or more samples missing"
fi
echo "**** Job ends ****"; date
exit $MISSING

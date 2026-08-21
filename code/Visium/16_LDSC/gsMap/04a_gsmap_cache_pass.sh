#!/bin/bash
#SBATCH --job-name=gsmap_cache
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cache_%A_%a.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cache_%A_%a.err
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-7
# -----------------------------------------------------------------------------
# 04a_gsmap_cache_pass.sh -- STAGE 1 of 2 for the full sweep.
#
# Builds the per-sample intermediates that every trait re-uses:
#   step 1 find_latent_representations (GVAE)   ~13 min
#   step 2 latent_to_gene (gene specificity)    ~19 min
#   step 3 generate_ldscore -> no-op in quick_mode (pre-built SNP x gene matrix)
#
# It runs ONE cheap trait (Height) purely to drive the pipeline through steps 1-3.
# Height is also the sweep's negative control, so this is not wasted work.
#
# After this finishes, 04b_gsmap_trait_array.sh fans out over the remaining traits
# and only pays step 4 (spatial_ldsc), reading these cached intermediates.
#
# Uses --gM_slices so gene ranks are on a common scale across sections
# => run 03_create_slice_mean.sh BEFORE this script.
#
#   sbatch 04a_gsmap_cache_pass.sh
# -----------------------------------------------------------------------------
set -uo pipefail

PROJ=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala
CODE=$PROJ/code/Visium/16_LDSC/gsMap
WORKDIR=$PROJ/processed-data/Visium/16_LDSC/gsMap
SCRATCH=/users/mtotty/claude_scratch/gsmap

ANNOT=BS_k16_Semisupervised_wAI
SLICE_MEAN=$WORKDIR/sample_slice_mean.parquet

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CODE/samples.txt)
mkdir -p "$CODE/logs"

echo "**** Job starts ****"; date
echo "sample: $SAMPLE   host: $(hostname)   cpus: ${SLURM_CPUS_PER_TASK:-NA}"

module load conda/3-24.3.0
source activate $SCRATCH/envs/gsmap

GM_ARG=""
if [ -f "$SLICE_MEAN" ]; then
    GM_ARG="--gM_slices $SLICE_MEAN"
    echo "using slice mean: $SLICE_MEAN"
else
    echo "WARNING: $SLICE_MEAN not found -- running WITHOUT --gM_slices."
    echo "         Gene ranks will not be comparable across sections."
fi

gsmap quick_mode \
    --workdir "$WORKDIR" \
    --sample_name "$SAMPLE" \
    --gsMap_resource_dir "$SCRATCH/gsMap_resource" \
    --hdf5_path "$WORKDIR/ST/${SAMPLE}.h5ad" \
    --annotation "$ANNOT" \
    --data_layer count \
    --sumstats_file "$WORKDIR/GWAS/height_ldscore.gz" \
    --trait_name Height \
    $GM_ARG \
    --max_processes ${SLURM_CPUS_PER_TASK:-8}

echo "**** Job ends ****"; date

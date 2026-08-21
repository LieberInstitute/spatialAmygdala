#!/bin/bash
#SBATCH --job-name=gsmap_sweep
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/gsmap_array_%A_%a.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/gsmap_array_%A_%a.err
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-7%3
# -----------------------------------------------------------------------------
# 04_gsmap_array.sh -- FULL SWEEP: 7 capture areas x 40 GWAS traits.
#
# One array task per capture area. Each task runs the complete gsMap pipeline
# (quick_mode = find_latent_representations -> latent_to_gene -> generate_ldscore
# -> spatial_ldsc -> cauchy_combination -> report) for ONE sample against ALL 40
# traits listed in gwas_config_full.yaml. The expensive per-sample stages (GVAE
# latent representation, GSS, LD scores) are computed once and reused across all
# 40 traits within the task, which is why traits are looped INSIDE a task rather
# than being their own array dimension.
#
# *** PREREQUISITES -- BOTH MUST COMPLETE FIRST ***
#   1. 03_create_slice_mean.sh   -> $WORKDIR/sample_slice_mean.parquet
#      Without it, --gM_slices has nothing to point at and gene ranks would be
#      section-local, making cross-sample results incomparable.
#   2. 03b_stage_gwas_symlinks.sh -> the 40 symlinks in $WORKDIR/GWAS/
#      Without them every path in gwas_config_full.yaml is dangling.
#   The script hard-fails below if either is missing.
#
# THROTTLING: --array=1-7%3 runs at most 3 tasks concurrently. Each task asks for
# 100G and 8 cpus, so an unthrottled 1-7 would hold 700G / 56 cores on `shared`
# at once -- that queues badly and is unneighbourly on a shared partition. %3
# caps the footprint at 300G / 24 cores while still finishing in roughly three
# waves. Drop to %2 if the partition is busy, or remove the %N if it is idle.
#
# MEMORY: gsMap docs quote ~80G for quick_mode at 120K cells; our sections are
# 26K-37K spots (Br6471 36876 ... Br9280 25835), so 100G is comfortably safe and
# leaves headroom for the 40-trait spatial_ldsc loop.
#
# --max_processes is set from $SLURM_CPUS_PER_TASK so gsMap's internal pool never
# oversubscribes the cpuset SLURM gave the task.
#
# Species is human -> NO --homolog_file anywhere.
# --data_layer count (SINGULAR): that is the layer name in our h5ad files; gsMap's
# default is "counts" and would fail.
#
# SUBMIT WITH (manual sbatch -- array jobs are not dispatched programmatically):
#   sbatch /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/04_gsmap_array.sh
#
# Re-run a single failed capture area, e.g. the 5th line of samples.txt:
#   sbatch --array=5 /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/04_gsmap_array.sh
# -----------------------------------------------------------------------------
set -uo pipefail

PROJ=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala
CODE=$PROJ/code/Visium/16_LDSC/gsMap
WORKDIR=$PROJ/processed-data/Visium/16_LDSC/gsMap
SCRATCH=/users/mtotty/claude_scratch/gsmap
RESOURCE=$SCRATCH/gsMap_resource

SAMPLE_FILE=$CODE/samples.txt
ANNOT=BS_k16_Semisupervised_wAI        # same domain labels as the classic LDSC run
DATA_LAYER=count                       # counts live in adata.layers['count'] (singular)
GM_SLICES=$WORKDIR/sample_slice_mean.parquet
GWAS_CONFIG=$CODE/gwas_config_full.yaml

mkdir -p "$CODE/logs"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_FILE")

echo "**** Job starts ****"; date
echo "host: $(hostname)  array job: ${SLURM_ARRAY_JOB_ID:-NA} task: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "sample: $SAMPLE  cpus: ${SLURM_CPUS_PER_TASK:-NA}  mem: ${SLURM_MEM_PER_NODE:-NA}"

# ---- preflight -------------------------------------------------------------
if [ -z "$SAMPLE" ]; then
    echo "FATAL: no sample on line ${SLURM_ARRAY_TASK_ID} of $SAMPLE_FILE"; exit 1
fi
if [ ! -f "$WORKDIR/ST/${SAMPLE}.h5ad" ]; then
    echo "FATAL: missing h5ad $WORKDIR/ST/${SAMPLE}.h5ad"; exit 1
fi
if [ ! -f "$GM_SLICES" ]; then
    echo "FATAL: missing $GM_SLICES -- run 03_create_slice_mean.sh first"; exit 1
fi
if [ ! -f "$GWAS_CONFIG" ]; then
    echo "FATAL: missing $GWAS_CONFIG"; exit 1
fi
# ----------------------------------------------------------------------------

module load conda/3-24.3.0
source activate $SCRATCH/envs/gsmap
echo "gsmap: $(which gsmap)  version: $(gsmap --version 2>&1 | head -1)"

gsmap quick_mode \
    --workdir "$WORKDIR" \
    --sample_name "$SAMPLE" \
    --gsMap_resource_dir "$RESOURCE" \
    --hdf5_path "$WORKDIR/ST/${SAMPLE}.h5ad" \
    --annotation "$ANNOT" \
    --data_layer "$DATA_LAYER" \
    --gM_slices "$GM_SLICES" \
    --sumstats_config_file "$GWAS_CONFIG" \
    --max_processes "${SLURM_CPUS_PER_TASK:-8}"

status=$?
echo "gsmap exit status: $status"

echo "---- outputs for $SAMPLE ----"
ls -la "$WORKDIR/$SAMPLE" 2>/dev/null | head -40

echo "**** Job ends ****"; date
exit $status

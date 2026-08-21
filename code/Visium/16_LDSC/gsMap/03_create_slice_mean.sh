#!/bin/bash
#SBATCH --job-name=gsmap_slice_mean
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/gsmap_slice_mean_%j.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/gsmap_slice_mean_%j.err
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
# -----------------------------------------------------------------------------
# 03_create_slice_mean.sh -- build the cross-section gene "slice mean" for gsMap.
#
# WHY THIS STEP EXISTS
#   gsMap's gene specificity score (GSS) is computed from a RANK of each gene's
#   expression within a spot's local neighbourhood. Computed independently per
#   capture area, those ranks are anchored to that one section's own expression
#   distribution: a gene that is mid-range in a section with deep sequencing can
#   land in a completely different rank bucket in a shallower section, purely
#   from a technical difference. The GSS values -- and therefore the spatial
#   LDSC chi^2 regressions built on top of them -- would then not be on a common
#   scale, and a domain-level p-value from Br6471 would not be comparable to the
#   same domain in Br9280.
#
#   `gsmap create_slice_mean` pools all seven sections and writes a single
#   per-gene mean expression reference (sample_slice_mean.parquet). Passing that
#   file to every quick_mode run via --gM_slices makes gsMap rank each gene
#   against the SHARED across-section reference instead of its own section, so
#   GSS -- and the LDSC results downstream -- are directly comparable across
#   capture areas, and the across-sample Cauchy combination in
#   05_cauchy_across_samples.sh is combining like with like.
#
#   This mirrors what the classic pipeline does by pseudobulking all sections
#   together before computing specificity (../LDSC/01_aggregate.R).
#
# INPUTS  : $WORKDIR/ST/<sample>.h5ad  (all 7 capture areas)
# OUTPUT  : $WORKDIR/sample_slice_mean.parquet
# RUN ORDER: this must finish BEFORE 04_gsmap_array.sh is submitted.
#
# NOTE on flags (verified against gsMap 1.73.7 --help, NOT the online docs):
#   create_slice_mean has NO --workdir; the output path is given in full by
#   --slice_mean_output_file. --data_layer is REQUIRED and must be `count`
#   (singular) because that is the layer name in our h5ad files -- gsMap's
#   default is `counts` (plural) and would fail.
#   Species is human, so NO --homolog_file.
#
# SUBMIT WITH:
#   sbatch /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/03_create_slice_mean.sh
# -----------------------------------------------------------------------------
set -uo pipefail

PROJ=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala
CODE=$PROJ/code/Visium/16_LDSC/gsMap
WORKDIR=$PROJ/processed-data/Visium/16_LDSC/gsMap
SCRATCH=/users/mtotty/claude_scratch/gsmap
RESOURCE=$SCRATCH/gsMap_resource

DATA_LAYER=count          # counts live in adata.layers['count'] (singular)

mkdir -p "$CODE/logs"
mkdir -p "$WORKDIR"

echo "**** Job starts ****"; date
echo "host: $(hostname)  cpus: ${SLURM_CPUS_PER_TASK:-NA}  mem: ${SLURM_MEM_PER_NODE:-NA}"

module load conda/3-24.3.0
source activate $SCRATCH/envs/gsmap
echo "gsmap: $(which gsmap)  version: $(gsmap --version 2>&1 | head -1)"

gsmap create_slice_mean \
    --sample_name_list \
        Br6471 \
        Br6660 \
        Br6423 \
        Br2743 \
        Br8325 \
        Br9192 \
        Br9280 \
    --h5ad_list \
        "$WORKDIR/ST/Br6471.h5ad" \
        "$WORKDIR/ST/Br6660.h5ad" \
        "$WORKDIR/ST/Br6423.h5ad" \
        "$WORKDIR/ST/Br2743.h5ad" \
        "$WORKDIR/ST/Br8325.h5ad" \
        "$WORKDIR/ST/Br9192.h5ad" \
        "$WORKDIR/ST/Br9280.h5ad" \
    --slice_mean_output_file "$WORKDIR/sample_slice_mean.parquet" \
    --data_layer "$DATA_LAYER"

echo "---- output ----"
ls -la "$WORKDIR/sample_slice_mean.parquet"

echo "**** Job ends ****"; date

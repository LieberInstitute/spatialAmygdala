#!/bin/bash
#SBATCH --job-name=gsmap_pilot
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/gsmap_pilot_%j.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/gsmap_pilot_%j.err
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=3-00:00:00
# -----------------------------------------------------------------------------
# 02_gsmap_quick_mode.sh -- gsMap PILOT: one capture area x 3 GWAS traits.
#
# Runs the whole gsMap pipeline in one command (quick_mode):
#   find_latent_representations (GVAE) -> latent_to_gene (GSS)
#   -> generate_ldscore -> spatial_ldsc -> cauchy_combination -> report
#
# quick_mode uses the pre-built 1000G EUR Phase3 SNP-by-gene weight matrix shipped
# in gsMap_resource/quick_mode/, so it needs ONLY --gsMap_resource_dir; the bfile,
# hapmap3 SNP list, GTF and LDSC weights are derived internally (gsMap/config.py).
#
# Genome build: the bundle GTF is gencode.v46lift37 = hg19, matching the hg19
# gene coordinates used by the classic pipeline (../LDSC/03_bed.R -> gene_meta_hg19.txt).
#
# Species is human, so NO --homolog_file.
#
# Intermediate outputs (latent representation, GSS, LD scores) are cached per sample
# in the workdir and are REUSED for any additional trait analysed later on this sample.
# -----------------------------------------------------------------------------
set -euo pipefail

PROJ=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala
CODE=$PROJ/code/Visium/16_LDSC/gsMap
WORKDIR=$PROJ/processed-data/Visium/16_LDSC/gsMap
SCRATCH=/users/mtotty/claude_scratch/gsmap
RESOURCE=$SCRATCH/gsMap_resource

SAMPLE=Br6471                            # pilot capture area (most spots, 14/15 domains)
ANNOT=BS_k16_Semisupervised_wAI           # same domain labels as the classic LDSC run
DATA_LAYER=count                          # counts live in layers['count'] (singular)

mkdir -p "$CODE/logs"

echo "**** Job starts ****"; date
echo "host: $(hostname)  cpus: ${SLURM_CPUS_PER_TASK:-NA}  mem: ${SLURM_MEM_PER_NODE:-NA}"

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
    --sumstats_config_file "$CODE/gwas_config_pilot.yaml" \
    --max_processes 8

echo "**** Job ends ****"; date

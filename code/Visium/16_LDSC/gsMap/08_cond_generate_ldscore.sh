#!/bin/bash
#SBATCH --job-name=gsmap_cond_ldsc
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_ldsc_%A_%a.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_ldsc_%A_%a.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-22
set -euo pipefail

# ARM is passed via --export=ARM=functional|neuronal
ARM="${ARM:?set ARM=functional or ARM=neuronal}"
SAMPLE="${SAMPLE:-Br6471}"
CHROM=$SLURM_ARRAY_TASK_ID

case "$ARM" in
  functional) ANNOT="/users/mtotty/claude_scratch/gsmap/additional_annotation/gsMap_additional_annotation" ;;
  neuronal)   ANNOT="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond/annotation_neuronal" ;;
  *) echo "bad ARM"; exit 2 ;;
esac

WORKDIR="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}"
mkdir -p "$WORKDIR"
export PATH="/users/mtotty/claude_scratch/gsmap/gsMap_resource/../envs/gsmap/bin:$PATH"

echo "ARM=$ARM SAMPLE=$SAMPLE CHROM=$CHROM WORKDIR=$WORKDIR"
echo "ANNOT=$ANNOT"
ls "$ANNOT/baseline.${CHROM}.annot.gz" || { echo "MISSING ANNOT"; exit 3; }

gsmap run_generate_ldscore \
  --workdir "$WORKDIR" \
  --sample_name "$SAMPLE" \
  --chrom "$CHROM" \
  --bfile_root "/users/mtotty/claude_scratch/gsmap/gsMap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC" \
  --keep_snp_root "/users/mtotty/claude_scratch/gsmap/gsMap_resource/LDSC_resource/hapmap3_snps/hm" \
  --gtf_annotation_file "/users/mtotty/claude_scratch/gsmap/gsMap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf" \
  --gene_window_size 50000 \
  --additional_baseline_annotation "$ANNOT"

echo "EXIT_OK chrom=$CHROM"
ls -la "$WORKDIR/$SAMPLE/generate_ldscore/" | head

#!/bin/bash
#SBATCH --job-name=neuro_annot
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond/logs/neuro_annot_%j.log
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond/logs/neuro_annot_%j.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --exclude=compute-058,compute-175
set -eo pipefail
echo "host=$(hostname) start=$(date)"
/users/mtotty/claude_scratch/gsmap/envs/gsmap/bin/python /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/build_annotation.py
echo "=== independent byte-identity verification ==="
OUT=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond/annotation_neuronal
for c in $(seq 1 22); do
  SRC=/users/mtotty/claude_scratch/gsmap/additional_annotation/gsMap_additional_annotation/baseline.$c.annot.gz
  NEW=$OUT/baseline.$c.annot.gz
  NC=$(zcat $NEW | head -1 | tr '\t' '\n' | wc -l)
  SC=$(zcat $SRC | head -1 | tr '\t' '\n' | wc -l)
  # strip the last (appended) field from every line of NEW and diff against SRC
  if zcat $NEW | cut -f1-$SC | cmp -s - <(zcat $SRC); then
     echo "chr$c VERIFY_OK src_cols=$SC new_cols=$NC $(zcat $NEW | wc -l) lines"
  else
     echo "chr$c VERIFY_FAIL"; exit 1
  fi
done
echo "end=$(date)"

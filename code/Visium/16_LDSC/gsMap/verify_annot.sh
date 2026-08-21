#!/bin/bash
#SBATCH --job-name=neuro_verify
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond/logs/neuro_verify_%j.log
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond/logs/neuro_verify_%j.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --exclude=compute-058,compute-175
set -e
OUT=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond/annotation_neuronal
SRCD=/users/mtotty/claude_scratch/gsmap/additional_annotation/gsMap_additional_annotation
echo "host=$(hostname) start=$(date)"
for c in $(seq 1 22); do
  SRC=$SRCD/baseline.$c.annot.gz
  NEW=$OUT/baseline.$c.annot.gz
  SC=$(zcat "$SRC" | head -1 | awk -F'\t' '{print NF}')
  NC=$(zcat "$NEW" | head -1 | awk -F'\t' '{print NF}')
  LASTNAME=$(zcat "$NEW" | head -1 | awk -F'\t' '{print $NF}')
  A=$(zcat "$SRC" | md5sum | awk '{print $1}')
  B=$(zcat "$NEW" | cut -f1-"$SC" | md5sum | awk '{print $1}')
  NL_S=$(zcat "$SRC" | wc -l)
  NL_N=$(zcat "$NEW" | wc -l)
  if [ "$A" == "$B" ] && [ "$NC" -eq "$((SC+1))" ] && [ "$LASTNAME" == "neuronal_axis" ] && [ "$NL_S" -eq "$NL_N" ]; then
     echo "chr$c VERIFY_OK src_cols=$SC new_cols=$NC lines=$NL_N md5=$A last_col=$LASTNAME"
  else
     echo "chr$c VERIFY_FAIL src_cols=$SC new_cols=$NC md5_src=$A md5_new=$B lines_src=$NL_S lines_new=$NL_N last=$LASTNAME"
  fi
done
echo "end=$(date)"

#!/bin/bash
#SBATCH --job-name=ctrl_b64
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/ctrl_b64.%j.log
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0:10:00
for f in fig_ctrl_heterogeneity.png fig_ctrl_donor_consistency.png fig_ctrl_confounds.png; do
  echo "@@@$f"
  base64 -w0 /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/scDRS/controls/$f
  echo
done

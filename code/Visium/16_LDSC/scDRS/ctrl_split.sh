#!/bin/bash
#SBATCH --job-name=ctrl_split
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/ctrl_split.%j.log
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0:15:00
mkdir -p /users/mtotty/claude_scratch/ctrl_xfer
cd /users/mtotty/claude_scratch/ctrl_xfer
for f in fig_ctrl_heterogeneity fig_ctrl_donor_consistency fig_ctrl_confounds; do
  base64 -w0 /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/scDRS/controls/$f.png > $f.b64
  split -b 50000 -d -a 3 $f.b64 $f.part
  echo "$f: $(stat -c %s $f.b64) b64 bytes, $(ls $f.part* | wc -l) parts"
done
ls -la /users/mtotty/claude_scratch/ctrl_xfer | head -5

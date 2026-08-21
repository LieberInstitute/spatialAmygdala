#!/bin/bash
#SBATCH --job-name=ctrl_split2
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/ctrl_split2.%j.log
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0:15:00
#
# Transfer helper. NEVER deletes: each invocation writes into its own
# subdirectory keyed by SLURM job id, so prior chunks are left in place.
set -euo pipefail
DEST="/users/mtotty/claude_scratch/ctrl_xfer/job_${SLURM_JOB_ID}"
mkdir -p "$DEST"
cd "$DEST"
for f in fig_ctrl_heterogeneity fig_ctrl_donor_consistency fig_ctrl_confounds; do
  base64 -w0 /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/scDRS/controls/$f.png > $f.b64
  split -b 45000 -d -a 3 $f.b64 $f.part
  echo "$f $(ls $f.part* | wc -l)"
done
echo "DEST=$DEST"

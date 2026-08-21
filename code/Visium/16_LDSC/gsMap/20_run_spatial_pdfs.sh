#!/bin/bash
#SBATCH --job-name=gsmap_spatial_pdfs
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --exclude=compute-053,compute-058,compute-175
#SBATCH --output=logs/spatial_pdfs_%j.out
#SBATCH --error=logs/spatial_pdfs_%j.err

set -eo pipefail
cd /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
GENV=/users/mtotty/claude_scratch/gsmap/envs/gsmap
echo "host $(hostname)  start $(date)"
$GENV/bin/python 20_spatial_pdfs.py
echo "end $(date)"
ls -la figures/spatial/ | tail -5

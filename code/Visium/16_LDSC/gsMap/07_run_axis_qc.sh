#!/bin/bash
#SBATCH --job-name=gsmap_axis_qc
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/axis_qc_%j.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/axis_qc_%j.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -eo pipefail
cd /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
echo "=== STAGE 1: pseudobulk ==="
/users/mtotty/claude_scratch/gsmap/gsMap_resource/../envs/gsmap/bin/python 07a_pseudobulk.py
echo; echo "=== STAGE 2: edgeR ==="
set +u; module load conda_R/4.5.x; set -u
Rscript 07b_neuronal_axis.R
echo "ALL DONE"

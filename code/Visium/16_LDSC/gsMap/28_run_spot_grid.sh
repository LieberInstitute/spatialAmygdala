#!/bin/bash
#SBATCH --job-name=spotgrid
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --exclude=compute-053,compute-058,compute-175
#SBATCH --output=logs/spotgrid_%j.out
#SBATCH --error=logs/spotgrid_%j.err
set -eo pipefail
cd /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
set +u; module load conda_R/4.5.x; set -u
Rscript 28_fig_spot_grid.R
echo "EXIT_OK"

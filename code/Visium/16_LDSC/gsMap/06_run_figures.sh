#!/bin/bash
#SBATCH --job-name=gsmap_figs
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/figs_%j.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/figs_%j.err
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
set -uo pipefail
CODE=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
module load conda/3-24.3.0
source activate /users/mtotty/claude_scratch/gsmap/envs/gsmap
python $CODE/06_pilot_figures.py
echo "**** done ****"

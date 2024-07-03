#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=02_build_spe
#SBATCH -o /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/visium_stitching/NacUtils/logs/02_build_spe_%a.log
#SBATCH -e /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/visium_stitching/NacUtils/logs/02_build_spe_%a.log
#SBATCH --array=6-6%5

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

line=$(awk -v task_id="$SLURM_ARRAY_TASK_ID" 'NR == task_id {print; exit}' sample_info.txt)

# Split the line into group and capture_area
group=$(echo "$line" | awk '{print $1}')

# Call the R script with the arrays
module load conda_R/devel

## List current modules for reproducibility
module list

Rscript 02_build_spe.R "$group"

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/
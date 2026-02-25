#!/bin/bash
#SBATCH --job-name=NMF_CV
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=32G
#SBATCH --output=logs/cv.txt
#SBATCH --error=logs/cv.txt
#SBATCH --mail-type=END
#SBATCH --time=72:00:00
#SBATCH --mail-user=mtotty2@jh.edu


echo "**** Job starts ****"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript NMF_cross_validation.R

echo "**** Job ends ****"
date
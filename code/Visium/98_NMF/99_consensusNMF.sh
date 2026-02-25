#!/bin/bash
#SBATCH --job-name=cNMF
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=150G
#SBATCH --output=logs/cNMF.txt
#SBATCH --error=logs/cNMF.txt
#SBATCH --mail-type=END
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
Rscript 99_consensusNMF.R

echo "**** Job ends ****"
date
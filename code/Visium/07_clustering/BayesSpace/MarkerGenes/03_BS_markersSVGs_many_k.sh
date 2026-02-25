#!/bin/bash
#SBATCH --job-name=BS_markersSVGs
#SBATCH --output=logs/BS_markersSVGs.txt
#SBATCH --error=logs/BS_markersSVGs.txt
#SBATCH --array=10-20
#SBATCH --mem=200G
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1 # specify the number of CPUs needed for the job, adjust as needed

echo "**** Job starts ****"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


## Load the R module
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 03_BS_markersSVGs_many_k.R

echo "**** Job ends ****"
date
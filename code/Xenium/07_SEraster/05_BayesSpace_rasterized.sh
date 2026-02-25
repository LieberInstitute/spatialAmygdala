#!/bin/bash
#SBATCH --job-name=BS_rasteri
#SBATCH --output=logs/bayesspace_%x.%a.txt
#SBATCH --error=logs/bayesspace_%x.%a.txt
#SBATCH --array=10-20
#SBATCH --mem=100G
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
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 05_BayesSpace_rasterized.R

echo "**** Job ends ****"
date
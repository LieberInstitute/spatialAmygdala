#!/bin/bash
#SBATCH --job-name=precast_many_k
#SBATCH --output=logs_precast/R-%x.%a.txt
#SBATCH --error=logs_precast/R-%x.%a.txt
#SBATCH --array=2-20
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1 # specify number of CPUs needed for the job, adjust as needed


echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript precast_many_k.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
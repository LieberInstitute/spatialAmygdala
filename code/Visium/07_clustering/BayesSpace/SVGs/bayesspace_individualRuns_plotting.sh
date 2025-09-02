#!/bin/bash
#SBATCH --job-name=individ_plotting
#SBATCH --output=logs/individualRuns_plotting.txt
#SBATCH --error=logs/individualRuns_plotting.txt
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
Rscript bayesspace_individualRuns_plotting.R

echo "**** Job ends ****"
date

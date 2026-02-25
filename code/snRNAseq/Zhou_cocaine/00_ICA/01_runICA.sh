#!/bin/bash
#SBATCH --job-name=ICA   # Job name
#SBATCH --output=./logs/ICA_Zhou.out  # Output file
#SBATCH --error=./logs/ICA_Zhou.err   # Error file
#SBATCH --ntasks=1                           
#SBATCH --mem=100G       
#SBATCH --cpus-per-task=30         # Number of CPU cores per task
#SBATCH --mail-type=END                         

echo "**** Job starts ****"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.5

## List current modules for reproducibility
module list

# Run the R script
Rscript 01_runICA.R

echo "**** Job ends ****"
date
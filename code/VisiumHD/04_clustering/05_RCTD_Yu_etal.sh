#!/bin/bash
#SBATCH --job-name=RCTD_HD   # Job name
#SBATCH --output=./logs/RCTD.out  # Output file
#SBATCH --error=./logs/RCTD.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=200G       
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
Rscript 05_RCTD_Yu_etal.R

echo "**** Job ends ****"
date
#!/bin/bash
#SBATCH --job-name=RCTD_Xenium_30   # Job name
#SBATCH --output=./logs/RCTD_30c.out  # Output file
#SBATCH --error=./logs/RCTD.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=300G       
#SBATCH --cpus-per-task=60         # Number of CPU cores per task
#SBATCH --mail-type=END            
#SBATCH --time=72:00:00             

echo "**** Job starts ****"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

## List current modules for reproducibility
module list

# Run the R script
Rscript 01_RCTD_updated.R

echo "**** Job ends ****"
date
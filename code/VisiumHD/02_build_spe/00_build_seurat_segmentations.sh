#!/bin/bash
#SBATCH --job-name=AmyHD_build_seurat   # Job name
#SBATCH --output=./logs/AmyHD_build_seurat.out  # Output file
#SBATCH --error=./logs/AmyHD_build_seurat.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=60G       
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
conda activate r-seurat-beta

## List current modules for reproducibility
module list

# Run the R script
Rscript 00_build_seurat_segmentations.R

echo "**** Job ends ****"
date
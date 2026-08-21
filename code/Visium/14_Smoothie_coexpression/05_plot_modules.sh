#!/bin/bash
#SBATCH --job-name=plot_modules
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/plot_modules.out
#SBATCH --error=logs/plot_modules.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=michael.totty@libd.org
 
set -euo pipefail
 
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node(s): ${SLURM_NODELIST}"
echo "Node memory requested: ${SLURM_MEM_PER_NODE}"
echo "n Tasks: ${SLURM_NTASKS}"
 
module load conda
conda activate smoothie_env
 
python 05_plot_modules.py
 
echo "Done: $(date)"
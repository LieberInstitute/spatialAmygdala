#!/bin/bash
#SBATCH --job-name=RCTD_Xenium
#SBATCH --output=./logs/RCTD_%a.out
#SBATCH --error=./logs/RCTD_%a.err
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=30
#SBATCH --mail-type=END
#SBATCH --time=72:00:00
#SBATCH --array=1-4

echo "**** Job starts ****"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURMD_NODENAME}"
echo "Array task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.4
module list

Rscript 02_RCTD_array_Totty.R ${SLURM_ARRAY_TASK_ID}

echo "**** Job ends ****"
date
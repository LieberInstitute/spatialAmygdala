#!/bin/bash
#SBATCH --job-name=RCTD_HD_array
#SBATCH --output=./logs/RCTD_%A_%a.out
#SBATCH --error=./logs/RCTD_%A_%a.err
#SBATCH --array=1-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=150G
#SBATCH --mail-type=END

echo "**** Job starts ****"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Array job id: ${SLURM_ARRAY_JOB_ID}"
echo "Array task id: ${SLURM_ARRAY_TASK_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURMD_NODENAME}"

module load conda_R/4.4
module list

Rscript 05_RCTD_Yu_ITCsubtypes_array.R

echo "**** Job ends ****"
date
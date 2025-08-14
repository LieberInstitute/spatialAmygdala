#!/bin/bash
#SBATCH --job-name=BANKSY
#SBATCH --output=logs/BANKSY_kmeans_%A_%a.txt
#SBATCH --error=logs/BANKSY_kmeans_%A_%a.txt
#SBATCH --mem=200G
#SBATCH --array=2-20
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mtotty2@jh.edu

echo "**** Job starts for k=${SLURM_ARRAY_TASK_ID} ****"
date

echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

# Load R
module load conda_R/4.4

# Run R script with the current value of k
Rscript banksy_many_k.R ${SLURM_ARRAY_TASK_ID}

echo "**** Job ends for k=${SLURM_ARRAY_TASK_ID} ****"
date

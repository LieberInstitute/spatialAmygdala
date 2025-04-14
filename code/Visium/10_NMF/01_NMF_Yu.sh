#!/bin/bash
#SBATCH --mem=32G
#SBATCH --job-name=NMF_Yu
#SBATCH -o logs/NMF_Yu.txt
#SBATCH -e logs/NMF_Yu.txt

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"


# load R
module load conda_R

# run R script
Rscript 01_NMF_Yu.R

# end script
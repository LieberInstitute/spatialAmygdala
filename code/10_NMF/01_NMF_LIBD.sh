#!/bin/bash
#SBATCH --mem=32G
#SBATCH --job-name=NMF_LIBD
#SBATCH -o logs/NMF_LIBD.txt
#SBATCH -e logs/NMF_LIBD.txt

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
Rscript 01_NMF_LIBD.R

# end script
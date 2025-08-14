#!/bin/bash
#SBATCH --job-name=Proseg_CCA
#SBATCH --output=logs/Proseg_CCA.out
#SBATCH --error=logs/Proseg_CCA.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"


module load conda_R/4.3

# Run the R script with an argument
Rscript seurat_cca.R 

echo "**** Job ends ****"
date


#!/bin/bash
#SBATCH --job-name=spatialAMY_dimred_GLM-PCA
#SBATCH --output=logs/spatialAMY_dimred_GLM-PCA.txt
#SBATCH --error=logs/spatialAMY_dimred_GLM-PCA.txt
#SBATCH --mem=250G
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=mtotty2@jh.edu


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"


## Load the R module 
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript dim_reduction.R

echo "**** Job ends ****"
date

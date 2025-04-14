#!/bin/bash
#SBATCH --job-name=spatialAMY_dimred
#SBATCH --output=logs/spatialAMY_dimred.txt
#SBATCH --error=logs/spatialAMY_dimred.txt
#SBATCH --mem=128G
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


## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.3

## List current modules for reproducibility
module list

## Edit with your job command
Rscript dim_reduction_pca.R

echo "**** Job ends ****"
date

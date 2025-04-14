#!/bin/bash
#SBATCH --job-name=AmyHD_banksy_8um
#SBATCH --output=logs/AmyHD_banksy_8um.txt
#SBATCH --error=logs/AmyHD_banksy_8um.txt
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


## Load the R module 
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript banksy_8um.R

echo "**** Job ends ****"
date

#!/bin/bash
#SBATCH --job-name=anksy
#SBATCH --output=logs/banksy.txt
#SBATCH --error=logs/anksy.txt
#SBATCH --mem=200G
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
module load conda_R/4.4

## Edit with your job command
Rscript banksy_testing.R

echo "**** Job ends ****"
date

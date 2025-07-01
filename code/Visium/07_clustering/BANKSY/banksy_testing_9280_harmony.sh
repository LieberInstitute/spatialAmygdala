#!/bin/bash
#SBATCH --job-name=BANKSY_9280
#SBATCH --output=logs/BANKSY_9280.txt
#SBATCH --error=logs/BANKSY_9280.txt
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
Rscript banksy_testing_9280_harmony.R

echo "**** Job ends ****"
date

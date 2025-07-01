#!/bin/bash
#SBATCH --job-name=build_xenium_15um
#SBATCH --output=logs/build_xenium_15um.txt
#SBATCH --error=logs/build_xenium_15um.txt
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=mtotty2@jh.edu
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
Rscript build_spe_xenium_15um.R

echo "**** Job ends ****"
date

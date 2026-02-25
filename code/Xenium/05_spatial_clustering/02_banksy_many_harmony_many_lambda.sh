#!/bin/bash
#SBATCH --job-name=Xenium_banksy
#SBATCH --output=logs/Xenium_banksy.txt
#SBATCH --error=logs/Xenium_banksy.txt
#SBATCH --mem=200G
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --time=72:00:00
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
Rscript 02_banksy_many_harmony_many_lambda.R

echo "**** Job ends ****"
date

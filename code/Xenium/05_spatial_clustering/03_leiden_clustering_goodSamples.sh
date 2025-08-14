#!/bin/bash
#SBATCH --job-name=Leiden_good
#SBATCH --output=logs/Leiden_good.txt
#SBATCH --error=logs/Leiden_good.txt
#SBATCH --mem=100G
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
module load conda_R/4.3

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 03_leiden_clustering_goodSamples.R

echo "**** Job ends ****"
date

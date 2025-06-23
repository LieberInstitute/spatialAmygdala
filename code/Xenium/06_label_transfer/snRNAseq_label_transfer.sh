#!/bin/bash
#SBATCH --job-name=label_transfer
#SBATCH --output=logs/snRNAseq_label_transfer.txt
#SBATCH --error=logs/snRNAseq_label_transfer.txt
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
Rscript snRNAseq_label_transfer.R

echo "**** Job ends ****"
date
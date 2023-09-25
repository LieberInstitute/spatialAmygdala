#!/bin/bash
#SBATCH --mem=150G
#SBATCH --job-name=BayesSpace_clustering
#SBATCH -o logs/BayesSpace_clustering.txt
#SBATCH -e logs/BayesSpace_clustering.txt
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript BayesSpace_clustering.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/

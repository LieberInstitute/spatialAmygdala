#!/bin/bash
#SBATCH --mem=150G
#SBATCH --cpus-per-task=10
#SBATCH --job-name=nnSVG_parallel
#SBATCH -o logs/nnSVG_parallel.txt
#SBATCH -e logs/nnSVG_parallel.txt
echo "**** Job starts ****"
datelog

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript nnSVG.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
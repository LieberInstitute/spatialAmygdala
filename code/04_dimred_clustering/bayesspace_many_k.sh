#!/bin/bash
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=60G
#$ -N bayesSpace_k_many
#$ -o logs/bayesSpace_k.$TASK_ID.txt
#$ -e logs/bayesSpace_k.$TASK_ID.txt
#$ -m e
#$ -t 5-20
#$ -tc 16

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript bayesspace_many_k.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
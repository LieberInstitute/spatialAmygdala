#!/bin/bash
#$ -cwd
#$ -l mem_free=400G,h_vmem=400G,h_fsize=10G
#$ -N Mouse_VisiumHD
#$ -o logs/clustering.txt
#$ -e logs/clustering.txt
#$ -m e

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
Rscript Seurat_clustering.R

echo "**** Job ends ****"
date
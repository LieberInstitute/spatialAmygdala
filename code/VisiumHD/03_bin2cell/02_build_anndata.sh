#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=64G
#SBATCH --job-name=02_basic_adata
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/02_basic_adata.%a.txt
#SBATCH -e logs/02_basic_adata.%a.txt
#SBATCH --array=1

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load bin2cell/0.3.0

## List current modules for reproducibility
module list

python 02_build_anndata.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
### This code was modified from /dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/code/04_bin2cell/
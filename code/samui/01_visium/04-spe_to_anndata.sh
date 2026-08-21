#!/bin/bash
#SBATCH --mem=40G
#SBATCH --job-name=04-spe2anndata
#SBATCH -o logs/spe2anndata_%a.txt
#SBATCH -e logs/spe2anndata_%a.txt
#SBATCH --array=1



echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.4
Rscript 03-spe_to_anndata.R

echo "**** Job ends ****"
date


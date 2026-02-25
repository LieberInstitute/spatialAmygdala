
#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=3G
#SBATCH --job-name=01_sample_info
#SBATCH -c 1
#SBATCH -t 10:00
#SBATCH -o logs/01_sample_info.txt
#SBATCH -e logs/01_sample_info.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_sample_info.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
### This code was modified from /dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/code/04_bin2cell/

#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=32G
#SBATCH --job-name=03_segment_and_b2c
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/03_segment_and_b2c.%a.txt
#SBATCH -e logs/03_segment_and_b2c.%a.txt
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

python 03_segment_and_b2c.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

#!/bin/bash
#SBATCH --mem=120G
#SBATCH -n 32
#SBATCH --job-name=segment_MeA
#SBATCH -o logs/segment_MeA.txt

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load spaceranger
export PATH=/users/$(whoami)/software/spaceranger/spaceranger-4.0.1:$PATH

module list

spaceranger segment --id=8325_MeA --tissue-image=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/images/VisiumHD/Br8325_AMY_anterior_MeA.tif

echo "**** Job ends ****"
date
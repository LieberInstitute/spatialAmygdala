#!/bin/bash
#SBATCH --mem=100G
#SBATCH -o svs2tif.txt 
#SBATCH -e svs2tif.txt
#SBATCH --array=1
#SBATCH --partition gpu

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"s
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load MATLAB
module load matlab/R2023a

matlab -nodesktop -nosplash -nojvm -r "run('/users/mtippani/Reviews/svs2tif.m')"

echo "**** Job ends ****"
date




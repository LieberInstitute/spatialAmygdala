#!/bin/bash
#SBATCH --job-name=HD_banksy_%a
#SBATCH --output=logs/HD_banksy_res%a.txt
#SBATCH --error=logs/HD_banksy_res%a.txt
#SBATCH --mem=200G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mtotty2@jh.edu
#SBATCH --array=0-5
 
echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Hostname: $(hostname)"
 
module load conda_R/4.5
module list
 
# Define resolutions to sweep
RESOLUTIONS=(0.2 0.4 0.6 0.8 1.0 1.2)
RES=${RESOLUTIONS[$SLURM_ARRAY_TASK_ID]}
 
echo "Running resolution: ${RES}"
Rscript 01_banksy_many_harmony_many_lambda.R ${RES}
 
echo "**** Job ends ****"
date
 
#!/bin/bash
#SBATCH --mem=90G
#SBATCH --job-name=01-create_samui
#SBATCH -o logs/samui_%a.txt
#SBATCH -e logs/samui_%a.txt
#SBATCH --array=1-5


# %4

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

donor=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/samui/visiumHD/sample_ids-visium.txt)
echo "Processing sample ${donor}"
date


module load samui/1.0.0-next.45
python 02-create_samui.py $donor

echo "**** Job ends ****"
date


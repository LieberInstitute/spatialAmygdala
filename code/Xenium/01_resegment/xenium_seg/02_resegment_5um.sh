#!/bin/bash
#SBATCH --job-name=5um_reseg
#SBATCH --output=logs/Xenium_resegment_5um.%A_%a.out
#SBATCH --error=logs/Xenium_resegment_5um.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --array=1

echo "**** Job starts ****"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Job Name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

## Load CellRanger
module load xeniumranger/2.0.0

## List current modules for reproducibility
module list

## Locate sample
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" reseg_last_sample_only.txt)
echo "Processing sample ${SAMPLE}"
echo "${SAMPLE}"
date

## Run CellRanger
xeniumranger resegment --id=${SAMPLE} \
                       --xenium-bundle=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/Xenium/copied_data/${SAMPLE} \
                       --localcores=64 \
                       --localmem=128 \
                       --expansion-distance=5

echo "**** Job ends ****"
date

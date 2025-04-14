#!/bin/bash
#SBATCH --job-name=Xenium_resegment
#SBATCH --output=logs/Xenium_resegment.%A_%a.out
#SBATCH --error=logs/Xenium_resegment.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=1-4%4

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
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" 01_resegment.txt)
echo "Processing sample ${SAMPLE}"
echo "${SAMPLE}"
date

## Run CellRanger
xeniumranger resegment --id=${SAMPLE} \
                       --xenium-bundle=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/Xenium/copied_data/2024-04-02_Psomagen/${SAMPLE} \
                       --localcores=16 \
                       --localmem=128 \
                       --expansion-distance=15

echo "**** Job ends ****"
date

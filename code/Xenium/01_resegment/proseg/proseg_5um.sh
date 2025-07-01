#!/bin/bash
#SBATCH --job-name=Proseg
#SBATCH --output=logs/Proseg.%A_%a.out
#SBATCH --error=logs/Proseg.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --array=1-4%4

echo "**** Job starts ****"
start_time=$(date +%s)
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Job Name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

## Load Rust
module load rust
export PATH="$HOME/.cargo/bin:$PATH"

## List current modules for reproducibility
module list

## Locate sample
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" 01_resegment.txt)
echo "Processing sample ${SAMPLE}"

## Set file paths/outs
INPUT_DIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Xenium/01_resegment/5um_segmentations/${SAMPLE}/outs
TRANSCRIPTS=${INPUT_DIR}/transcripts.csv.gz

## Create output directory if needed
OUTPUT_DIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Xenium/01_resegment/proseg_5um/${SAMPLE}
mkdir -p ${OUTPUT_DIR}

## Run proseg
proseg --xenium --output-path=${OUTPUT_DIR} ${TRANSCRIPTS}


echo "**** Job ends ****"
end_time=$(date +%s)
date

## Calculate and print runtime
runtime=$((end_time - start_time))
echo "Total runtime: ${runtime} seconds"

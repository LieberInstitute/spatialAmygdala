#!/bin/bash
#SBATCH --job-name=SegCompare
#SBATCH --output=logs/SegCompare_%A_%a.out
#SBATCH --error=logs/SegCompare_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --array=1-4

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"


module load conda_R/4.3

# Define the list of R objects to process
OBJ_LIST=("spe_xenium_5um" "spe_xenium_15um" "spe_proseg_5um" "spe_proseg_15um")

# Select based on array task ID
OBJ_NAME=${OBJ_LIST[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing $OBJ_NAME"

# Run the R script with an argument
Rscript comparing_segmentations_good_samples.R $OBJ_NAME

echo "**** Job ends ****"
date


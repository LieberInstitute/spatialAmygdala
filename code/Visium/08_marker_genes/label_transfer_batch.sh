#!/bin/bash
#SBATCH --job-name=LabelTransfer
#SBATCH --output=logs/label_transfer_%A_%a.out
#SBATCH --error=logs/label_transfer_%A_%a.err
#SBATCH --array=1-9
#SBATCH --mem=64G

# Your sample list (excluding Br8325)
SAMPLES=(Br2743 Br6423 Br6471 Br6660 Br9017 Br9192 Br9206 Br9280 Br9469)

QUERY_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

echo "Running label transfer: Br8325 → $QUERY_ID"

# load R
module load conda_R/4.4

# Run your R script for this sample
Rscript label_transfer_batch.R $QUERY_ID

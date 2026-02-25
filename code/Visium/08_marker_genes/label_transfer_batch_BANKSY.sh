#!/bin/bash
#SBATCH --job-name=BANKSYTransfer
#SBATCH --mem=64G
#SBATCH --array=1-9
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/banksy_transfer_%A_%a.out
#SBATCH --error=logs/banksy_transfer_%A_%a.err

# Sample list (excluding Br8325)
SAMPLES=(Br2743 Br6423 Br6471 Br6660 Br9017 Br9192 Br9206 Br9280 Br9469)
QUERY_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# Static log files per sample


echo "Running label transfer: Br8325 → $QUERY_ID"

# load R
module load conda_R/4.4

# Run your R script for this sample
Rscript label_transfer_batch_BANKSY.R $QUERY_ID

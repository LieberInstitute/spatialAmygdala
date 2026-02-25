#!/bin/bash
#$ -cwd
#$ -l mem_free=400G
#$ -N subset_ann_data
#$ -o logs/subset_ann_data.out
#$ -e logs/subset_ann_data.err
#$ -m a

# Activate conda (optional)
module load conda_R/4.4
conda activate allen_env

# Run Python script
python loading_h5ad_AnnData.py
#!/bin/bash
#SBATCH --job-name=fix_coords
#SBATCH --output=logs/fix_banksy_coords_%j.out
#SBATCH --error=logs/fix_banksy_coords_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mtotty2@jh.edu

echo "**** Job starts ****"
date
echo "Running on: $HOSTNAME"
echo "User: $USER"

# Load R
module load conda_R/4.4
module list

# Run the coordinate reset script
Rscript 08_fix_coords.R

echo "**** Job ends ****"
date

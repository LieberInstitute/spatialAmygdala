#!/bin/bash
#SBATCH --job-name=plot_banksy
#SBATCH --output=logs/plot_banksy_%j.out
#SBATCH --error=logs/plot_banksy_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mtotty2@jh.edu

echo "**** Job starts ****"
date
echo "Running on: $HOSTNAME"
echo "User: $USER"

# Load R module
module load conda_R/4.4
module list

# Run the plotting script
Rscript 08_plot_batch_xenium.R

echo "**** Job ends ****"
date

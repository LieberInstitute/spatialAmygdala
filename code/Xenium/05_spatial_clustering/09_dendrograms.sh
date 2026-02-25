#!/bin/bash
#SBATCH --job-name=banksy_dendros
#SBATCH --output=logs/plot_banksy_dendros_%j.out
#SBATCH --error=logs/plot_banksy_dendros_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mtotty2@jh.edu

echo "**** Job starts ****"
date
module load conda_R/4.4
module list

Rscript 09_dendrograms.R

echo "**** Job ends ****"
date
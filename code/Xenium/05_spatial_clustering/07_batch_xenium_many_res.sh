#!/bin/bash
#SBATCH --job-name=Xenium_banksy_resSweep
#SBATCH --output=logs/Xenium_banksy_resSweep_%A_%a.out
#SBATCH --error=logs/Xenium_banksy_resSweep_%A_%a.err
#SBATCH --array=0-7
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mtotty2@jh.edu

echo "**** Job starts ****"
date

module load conda_R/4.4
module list

# Resolution sweep
RESOLUTIONS=(0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0)
RES=${RESOLUTIONS[$SLURM_ARRAY_TASK_ID]}

Rscript 07_batch_xenium_many_res.R $RES

echo "**** Job ends ****"
date

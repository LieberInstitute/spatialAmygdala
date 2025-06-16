#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=0-05:00:00
#SBATCH --array=1-10
#SBATCH --job-name=nnSVG_perSample
#SBATCH --output=logs/nnSVG_perSample_%A_%a.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mtotty2@jh.edu

echo "**** Job starts ****"
date

input=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 01_nnSVG_arrayjobs.txt)
echo "Running on sample: $input"

module load conda_R/4.4
module list

Rscript 01_nnSVG_arrayjobs.R $input

echo "**** Job ends ****"
date

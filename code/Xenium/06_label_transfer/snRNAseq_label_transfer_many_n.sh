#!/bin/bash
#SBATCH --job-name=label_transfer
#SBATCH --output=logs/label_transfer_%A_%a.out
#SBATCH --error=logs/label_transfer_%A_%a.err
#SBATCH --array=1-7
#SBATCH --mem=128G
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mtotty2@jh.edu

module load conda_R/4.4
module list

# Define array of de.n values
de_n_values=(3 5 10 15 25 35 50)

# Get de.n corresponding to SLURM_ARRAY_TASK_ID
de_n=${de_n_values[$SLURM_ARRAY_TASK_ID-1]}

echo "**** Running SingleR with de.n = ${de_n} ****"
Rscript snRNAseq_label_transfer_many_n.R $de_n
echo "**** Job complete for de.n = ${de_n} ****"

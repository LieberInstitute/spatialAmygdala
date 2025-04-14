#!/bin/bash
#SBATCH --job-name=banksy_lambda
#SBATCH --output=logs/banksy_lambda_%a.out
#SBATCH --error=logs/banksy_lambda_%a.err
#SBATCH --array=1-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=02:00:00

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.4 

# Define lambda values
LAMBDAS=(0.2 0.4 0.6 0.8)

# Get the lambda for this job array index
LAMBDA=${LAMBDAS[$SLURM_ARRAY_TASK_ID-1]}

# Print lambda for debugging
echo "Running with lambda=${LAMBDA}"

# Run the R script with the current lambda value
Rscript banksy_many_lambda.R $LAMBDA

echo "**** Job ends ****"
date


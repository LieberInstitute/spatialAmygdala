#!/bin/bash
#SBATCH --mem=80G
#SBATCH -n 8
#SBATCH --job-name=AMY-HD_spaceranger
#SBATCH --output=logs/HD-spaceranger-2509-%a.txt
#SBATCH --array=1-4


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load spaceranger
module load spaceranger/4.0.1


## Locate file
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' 03-sample_ids-2509.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
IMAGE=$(awk 'BEGIN {FS="\t"} {print $2}' 03-sample_ids-2509.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
IMGCYT=$(awk 'BEGIN {FS="\t"} {print $3}' 03-sample_ids-2509.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SAMPARG=$(awk 'BEGIN {FS="\t"} {print $4}' 03-sample_ids-2509.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
# SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" 03_24-09_samples-list.txt)
echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDE=$(echo ${SAMPLE} | cut -d "_" -f 1)
CAPTUREAREA=$(echo ${SAMPLE} | cut -d "_" -f 2)
SAM=$(paste <(echo ${SLIDE}) <(echo "-") <(echo ${CAPTUREAREA}) -d '')
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}"

## Find FASTQ file path
FASTQPATH=$(ls -d /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/FASTQ/VisiumHD/)

## Hank from 10x Genomics recommended setting this environment
export NUMBA_NUM_THREADS=1

spaceranger count \
    --id=${SAMPLE} \
    --sample=${SAMPARG} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATH} \
    --probe-set=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/HD_probe_set/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --cytaimage=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/images/VisiumHD/${IMGCYT}.tif \
    --image=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/images/VisiumHD/${IMAGE}.tif \
    --create-bam=false \
    --localcores=8 \
    --localmem=64 


## Move output
echo "Moving results to new location"
date
mkdir -p /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger
mv ${SAMPLE} /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger

echo "**** Job ends ****"
date

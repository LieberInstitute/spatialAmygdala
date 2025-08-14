#!/bin/bash
#SBATCH --mem=80G
#SBATCH -n 32
#SBATCH --job-name=HD_spaceranger
#SBATCH -o logs/HD-spaceranger-CEA%a.o.txt
#SBATCH --array=1

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load SpaceRanger
module load spaceranger/3.0.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" sample_ids.txt)
echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDE=$(echo ${SAMPLE} | cut -d "_" -f 1)
CAPTUREAREA=$(echo ${SAMPLE} | cut -d "_" -f 2)
SAM=$(paste <(echo ${SLIDE}) <(echo "-") <(echo ${CAPTUREAREA}) -d '')
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}"

## Find FASTQ file path
FASTQPATH=$(ls -d /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/FASTQ/VisiumHD/${SAMPLE}/)


## Hank from 10x Genomics recommended setting this environment
export NUMBA_NUM_THREADS=1

spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATH} \
    --probe-set=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/HD_probe_set/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --cytaimage=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/images/VisiumHD/CAVG10676_2024-08-21_13-19-15_2024-08-21_12-55-59_H1-W369TJK_A1_HDp_s003_AMY.tif \
    --image=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/images/VisiumHD/40x_HE_AMY_s003.tif \
    --create-bam=false \
    --localcores=32 \
    --localmem=64 \


## Move output
echo "Moving results to new location"
date
mkdir -p /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/01_spaceranger
mv ${SAMPLE} /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/01_spaceranger

echo "**** Job ends ****"
date
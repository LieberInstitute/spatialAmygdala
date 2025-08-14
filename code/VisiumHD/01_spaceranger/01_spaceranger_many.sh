#!/bin/bash
#SBATCH --mem=80G
#SBATCH --cpus-per-task=32
#SBATCH --job-name=HD_spaceranger
#SBATCH --output=logs/HD-spaceranger-%a.out
#SBATCH --error=logs/HD-spaceranger-%a.err
#SBATCH --array=1-5


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load spaceranger
export PATH=/users/$(whoami)/software/spaceranger/spaceranger-4.0.1:$PATH

module list

## Parse sample info
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" sample_ids.txt)
echo "Processing sample ${SAMPLE}"

SLIDE=$(echo ${SAMPLE} | cut -d "_" -f 1)
CAPTUREAREA=$(echo ${SAMPLE} | cut -d "_" -f 2)
echo "Slide: ${SLIDE}, Area: ${CAPTUREAREA}"

## Find FASTQ
FASTQPATH=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/FASTQ/VisiumHD/${SAMPLE}/

## Define image directory
IMAGEDIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/images/VisiumHD

## Dynamically find the cytassist image
CYTAIMAGE=$(find ${IMAGEDIR} -type f -iname "*${SAMPLE}*_HDp*.tif" | head -n 1)
if [[ -z "${CYTAIMAGE}" ]]; then
    echo "ERROR: Cytassist image not found for ${SAMPLE}"
    exit 1
fi

## Dynamically find the high-res HE image (exclude HDp to avoid matching cyta)
IMAGE=$(find ${IMAGEDIR} -type f -iname "${SAMPLE}*.tif" ! -iname "*HDp*" | head -n 1)
if [[ -z "${IMAGE}" ]]; then
    echo "ERROR: HE image not found for ${SAMPLE}"
    exit 1
fi

echo "Cytassist image: ${CYTAIMAGE}"
echo "HE image: ${IMAGE}"

## Set Numba environment
export NUMBA_NUM_THREADS=1

## Run spaceranger
spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATH} \
    --probe-set=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/raw-data/HD_probe_set/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --cytaimage=${CYTAIMAGE} \
    --image=${IMAGE} \
    --create-bam=false \
    --localcores=32 \
    --localmem=64

## Move output
echo "Moving results to ${OUTPUTDIR}"
date
OUTPUTDIR=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/01_spaceranger/${SAMPLE}
mkdir -p ${OUTPUTDIR}
mv ${SAMPLE} ${OUTPUTDIR}


echo "**** Job ends ****"
date

#!/bin/bash
#$ -cwd
#$ -pe local 10
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o logs/$TASK_ID_splitSlide.txt
#$ -e logs/$TASK_ID_splitSlide.txt
#$ -t 1

 
echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"

TOOLBOX='/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VistoSeg/code/'
IMAGELOC='/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/raw-data/Images/ifImages/V12N28-333' 
FNAME=${IMAGELOC}/V12N28-333.mat
SLIDEID='V12N28-333'

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$TOOLBOX')), splitSlide_IF('$FNAME')"

mkdir -p /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/Images/VistoSeg/Capture_areas/${SLIDEID}/

mv ${IMAGELOC}/${SLIDEID}_*1* /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/Images/VistoSeg/Capture_areas/${SLIDEID}/

echo "**** Job ends ****"
date


#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o logs/InFormStitch.txt
#$ -e logs/InFormStitch.txt

 
echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"

module load matlab/R2019a

toolbox='/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VistoSeg/code/'

filename='/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/raw-data/Images/ifImages/V12N28-333/20230524_VSPG_dACC_4celltypes_Scan1_*_component_data.tif'
fname='V12N28-333'
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O{1} = 'DAPI'; O{2} = 'NeuN'; O{3} = 'TMEM119'; O{4} = 'GFAP'; O{5} = 'OLIG2'; O{6} = 'AF'; InFormStitch('$filename',O,6,'$fname')"


echo "**** Job ends ****"
date

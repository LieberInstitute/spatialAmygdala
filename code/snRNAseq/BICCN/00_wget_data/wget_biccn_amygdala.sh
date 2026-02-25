#!/bin/bash
#SBATCH --mem=5G
#SBATCH --job-name=BICCN_AMY
#SBATCH -o logs/BICCN_AMY_download.txt
#SBATCH -e logs/BICCN_AMY_download.txt

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

# URL of the file to be downloaded
central_url="https://datasets.cellxgene.cziscience.com/9287f4c4-3764-4d33-a905-0a9bdbc064e2.rds"
medial_url="https://datasets.cellxgene.cziscience.com/bb88b65d-978a-47ba-9cfc-b227cc6b2516.rds"
basolateral_url="https://datasets.cellxgene.cziscience.com/4e124ecc-7885-465c-bab9-4e94d9d40b6a.rds"
basomedial_url="https://datasets.cellxgene.cziscience.com/5bf5d239-6699-4df5-a9cf-406f533fc178.rds"
lateral_url="https://datasets.cellxgene.cziscience.com/fa820927-923b-4a90-86a1-7ce0b6b4335f.rds"
cortical_url="https://datasets.cellxgene.cziscience.com/96c467f9-e0c5-4d30-98e9-57aefdfd5ddb.rds"


# Target file name to save as (optional, can be changed to your preference)
central_fname="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/snRNA-seq/BICCN/biccnAMY_central.rds"
medial_fname="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/snRNA-seq/BICCN/biccnAMY_medial.rds"
basolateral_fname="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/snRNA-seq/BICCN/biccnAMY_basolateral.rds"
basomedial_fname="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/snRNA-seq/BICCN/biccnAMY_basomedial.rds"
lateral_fname="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/snRNA-seq/BICCN/biccnAMY_lateral.rds"
cortical_fname="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/snRNA-seq/BICCN/biccnAMY_cortical.rds"



# Using wget to download the file
echo "Starting $central_fname download"
wget -O $central_fname $central_url

echo "Starting $medial_fname download"
wget -O $medial_fname $medial_url

echo "Starting $basolateral_fname download"
wget -O $basolateral_fname $basolateral_url
w
echo "Starting $basomedial_fname donload"
wget -O $basomedial_fname $basomedial_url

echo "Starting $lateral_fname download"
wget -O $lateral_fname $lateral_url

echo "Starting $cortical_fname download"
wget -O $cortical_fname $cortical_url

# End of script
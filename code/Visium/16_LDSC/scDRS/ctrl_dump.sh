#!/bin/bash
#SBATCH --job-name=ctrl_dump
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/ctrl_dump.%j.log
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0:10:00
echo '@@@ctrl_stats_summary.tsv'; cat /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/scDRS/controls/ctrl_stats_summary.tsv
echo '@@@ctrl_heterogeneity.tsv'; cat /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/scDRS/controls/ctrl_heterogeneity.tsv
echo '@@@ctrl_domain_covariates.tsv'; cat /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/scDRS/controls/ctrl_domain_covariates.tsv
echo '@@@ctrl_gwas_power.tsv'; cat /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/scDRS/controls/ctrl_gwas_power.tsv
echo '@@@ctrl_spot_covariates.tsv'; cat /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/scDRS/controls/ctrl_spot_covariates.tsv

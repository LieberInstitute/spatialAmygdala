#!/bin/bash
# -----------------------------------------------------------------------------
# 03b_stage_gwas_symlinks.sh -- stage the 40 classic-LDSC GWAS sumstats into
# $WORKDIR/GWAS/ as SYMLINKS.
#
# gsMap wants every sumstats file reachable from one directory named in
# gwas_config_full.yaml. The authoritative files are read-only shared data at
#   /dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain/
# so we symlink rather than copy: no duplicated multi-GB files, and the source
# is never touched.
#
# IDEMPOTENT: every link is created only if the path does not already exist,
# so re-running is harmless. Three links (SCZ / MDD / Height) already exist
# from the pilot and are left alone. NOTHING IS EVER DELETED OR OVERWRITTEN.
#
# This is a light file-op script -- run it on the login node, no sbatch needed:
#   bash /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/03b_stage_gwas_symlinks.sh
# -----------------------------------------------------------------------------
set -uo pipefail

PROJ=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala
CODE=$PROJ/code/Visium/16_LDSC/gsMap
WORKDIR=$PROJ/processed-data/Visium/16_LDSC/gsMap
SCRATCH=/users/mtotty/claude_scratch/gsmap
RESOURCE=$SCRATCH/gsMap_resource

GWAS_SRC=/dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain
GWAS_DST=$WORKDIR/GWAS

echo "**** Job starts ****"; date

mkdir -p "$GWAS_DST"

n_new=0
n_have=0
n_miss=0

link_one() {
    local f="$1"
    if [ ! -e "$GWAS_SRC/$f" ]; then
        echo "MISSING SOURCE: $GWAS_SRC/$f"
        n_miss=$((n_miss+1))
        return
    fi
    if [ -e "$GWAS_DST/$f" ]; then
        echo "exists : $f"
        n_have=$((n_have+1))
    else
        ln -s "$GWAS_SRC/$f" "$GWAS_DST/$f"
        echo "linked : $f"
        n_new=$((n_new+1))
    fi
}

link_one adhd.gz
link_one alcohol_ldscore.gz
link_one Alzheimer_ldscore.gz
link_one Alzheimer_ldscore2.gz
link_one Alzheimer_ldscore3.gz
link_one Anorexia_ldscore.gz
link_one AUD.EUR_META.NatMed2023.gz
link_one EA_aud_Jun07.txt_ldscore.gz
link_one autism_ldscore.gz
link_one bp_ldscore.gz
link_one bp3_ldscore.gz
link_one bmi_ldscore.gz
link_one UKB_460K.cov_EDU_YEARS.sumstats.gz
link_one epilepsyAll.gz
link_one epilepsyFocal.gz
link_one epilepsyGGE.gz
link_one GSCAN_AgeSmk_2022_ldscore.gz
link_one GSCAN_CigDay_2022_ldscore.gz
link_one GSCAN_DrnkWk_2022_ldscore.gz
link_one GSCAN_SmkCes_2022_ldscore.gz
link_one GSCAN_SmkInit_2022_ldscore.gz
link_one height_ldscore.gz
link_one insomnia.gz
link_one intelligence.gz
link_one MDD_ldscore_PGC_UKB_23andme.gz
link_one mdd2019edinburgh_ldscore.gz
link_one mdd_ex23andMe_ldscore.gz
link_one UKB_460K.mental_NEUROTICISM.sumstats.gz
link_one OUD_EA_MVP1_Mar12.gz
link_one OUD_EA_MVP2_Mar12.gz
link_one OUD_EA_MVP1_MVP2_YP_SAGE_Mar12.gz
link_one PD_ldscore.gz
link_one PTSD.gz
link_one scz_PGC3_ldscore.gz
link_one scz_PGC2_CLOZUK_ldscore.gz
link_one smoking_ldscore.gz
link_one MEGASTROKE.1.AS.EUR_ldscore.gz
link_one MEGASTROKE.2.AIS.EUR_ldscore.gz
link_one stroke2022any.gz
link_one PASS_Type_2_Diabetes.sumstats.gz

echo "----"
echo "newly linked : $n_new"
echo "already there: $n_have"
echo "missing src  : $n_miss"
echo "total links in $GWAS_DST: $(ls -1 "$GWAS_DST" | wc -l)"

echo "**** Job ends ****"; date

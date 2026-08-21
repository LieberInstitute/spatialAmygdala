#!/bin/bash
#SBATCH --job-name=gsmap_cauchy_all
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/gsmap_cauchy_all_%j.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/gsmap_cauchy_all_%j.err
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --time=3-00:00:00
# -----------------------------------------------------------------------------
# 05_cauchy_across_samples.sh -- domain-level p-values pooled ACROSS all 7 sections.
#
# quick_mode already runs cauchy_combination per sample, giving one p-value per
# (domain, trait, SAMPLE). This script instead pools the spot-level p-values from
# ALL SEVEN capture areas into ONE p-value per (domain, trait), which is the
# quantity that lines up with the classic 15-domain x 40-trait ldsc_results.csv
# (the classic run pseudobulks all sections together, so its unit is the domain,
# not the domain-within-section).
#
# The Cauchy combination test is used because spot-level p-values within a domain
# are strongly dependent (neighbouring spots share LD-score signal); ACAT/Cauchy
# is valid under arbitrary dependence, unlike Fisher's method.
#
# Pooling across samples is only meaningful because every quick_mode task in
# 04_gsmap_array.sh ranked genes against the SHARED slice mean
# ($WORKDIR/sample_slice_mean.parquet) -- see 03_create_slice_mean.sh.
#
# *** PREREQUISITE: 04_gsmap_array.sh must have completed for ALL 7 samples. ***
#
# Flags (verified against gsMap 1.73.7 --help):
#   run_cauchy_combination requires --workdir --trait_name --annotation;
#   --sample_name_list takes space-separated sample names and --output_file is
#   REQUIRED whenever more than one sample is combined.
#
# Trait names are read straight out of gwas_config_full.yaml so the two files can
# never drift apart.
#
# OUTPUT: $WORKDIR/cauchy_across_samples/<trait>_cauchy.csv.gz  (one per trait)
#
# SUBMIT WITH:
#   sbatch /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/05_cauchy_across_samples.sh
#
# Or chain it directly behind the array (replace <ARRAYJOBID>):
#   sbatch --dependency=afterok:<ARRAYJOBID> /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/05_cauchy_across_samples.sh
# -----------------------------------------------------------------------------
set -uo pipefail

PROJ=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala
CODE=$PROJ/code/Visium/16_LDSC/gsMap
WORKDIR=$PROJ/processed-data/Visium/16_LDSC/gsMap
SCRATCH=/users/mtotty/claude_scratch/gsmap
RESOURCE=$SCRATCH/gsMap_resource

ANNOT=BS_k16_Semisupervised_wAI
GWAS_CONFIG=$CODE/gwas_config_full.yaml
SAMPLE_FILE=$CODE/samples.txt
OUTDIR=$WORKDIR/cauchy_across_samples

mkdir -p "$CODE/logs"
mkdir -p "$OUTDIR"

echo "**** Job starts ****"; date
echo "host: $(hostname)  cpus: ${SLURM_CPUS_PER_TASK:-NA}  mem: ${SLURM_MEM_PER_NODE:-NA}"

module load conda/3-24.3.0
source activate $SCRATCH/envs/gsmap
echo "gsmap: $(which gsmap)  version: $(gsmap --version 2>&1 | head -1)"

# all 7 capture areas, space separated, from samples.txt
SAMPLE_LIST=$(tr '\n' ' ' < "$SAMPLE_FILE")
echo "samples: $SAMPLE_LIST"

# trait names = the yaml keys (skip comments and blank lines)
TRAITS=$(grep -v '^[[:space:]]*#' "$GWAS_CONFIG" | grep ':' | sed 's/:.*//' | sed 's/[[:space:]]*$//' | grep -v '^$')
echo "n traits: $(echo "$TRAITS" | wc -l)"

n_ok=0
n_fail=0

for TRAIT in $TRAITS; do
    OUT=$OUTDIR/${TRAIT}_cauchy.csv.gz
    echo "==== $TRAIT ===="; date
    gsmap run_cauchy_combination \
        --workdir "$WORKDIR" \
        --trait_name "$TRAIT" \
        --annotation "$ANNOT" \
        --sample_name_list $SAMPLE_LIST \
        --output_file "$OUT"
    if [ $? -eq 0 ] && [ -s "$OUT" ]; then
        echo "OK   $TRAIT -> $OUT"
        n_ok=$((n_ok+1))
    else
        echo "FAIL $TRAIT"
        n_fail=$((n_fail+1))
    fi
done

echo "----"
echo "succeeded: $n_ok   failed: $n_fail"
ls -la "$OUTDIR"

echo "**** Job ends ****"; date

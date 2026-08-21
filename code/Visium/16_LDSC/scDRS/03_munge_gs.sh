#!/bin/bash
#SBATCH --job-name=scdrs_munge
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/03_munge_gs.%j.log
# -----------------------------------------------------------------------------
# 03_munge_gs.sh -- turn the MAGMA gene-level results into scDRS gene sets.
#
#   sbatch 03_munge_gs.sh          (run AFTER 02_magma_array.sh finishes)
#
#   1. 03a_build_zscore_matrix.py  -> $GS_DIR/magma_zscore_matrix.tsv
#   2. scdrs munge-gs              -> $GS_DIR/magma_top${GS_NMAX}.gs
#
# munge-gs takes the top --n-max genes per trait by Z and writes one line per
# trait: TRAIT <tab> GENE:WEIGHT,GENE:WEIGHT,...  Weights are the Z-scores
# (--weight zscore), so a gene's contribution to the disease score scales with
# its GWAS evidence rather than every gene counting equally.
#
# n-max = 1000 matches the scDRS paper's published gene sets.
# -----------------------------------------------------------------------------
set -euo pipefail
# NOTE: under sbatch, SLURM copies this script to a spool directory, so "$0"
# does NOT resolve to this file's real location -- $(dirname $(readlink -f $0))
# yields /var/spool/slurm/... and config.sh is not found there. Hardcode the
# code directory instead. Override with CODE_DIR=... if you relocate the tree.
CODE_DIR=${CODE_DIR:-/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS}
source "$CODE_DIR/config.sh"

echo "**** Job starts ****"; date; hostname
activate_env
mkdir -p "$GS_DIR" "$CODE/logs"

ZMAT=$GS_DIR/magma_zscore_matrix.tsv
GS=$GS_DIR/magma_top${GS_NMAX}.gs

echo "[1/2] assembling gene x trait Z matrix"
"$PY" "$CODE/03a_build_zscore_matrix.py" "$MAGMA_OUT" "$TRAITS_TSV" "$ZMAT" \
    | tee "$CODE/03_zscore_matrix_summary.txt"

echo
echo "[2/2] scdrs munge-gs"
"$SCDRS" munge-gs \
    --out-file "$GS" \
    --zscore-file "$ZMAT" \
    --weight zscore \
    --n-max "$GS_NMAX"

echo "---- $GS ----"
cut -f1 "$GS" | head -50
echo "trait lines: $(( $(wc -l < "$GS") - 1 ))"

echo "**** Job ends ****"; date

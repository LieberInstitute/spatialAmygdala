#!/bin/bash
#SBATCH --job-name=scdrs_magma
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-40%20
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/02_magma.%A_%a.log
# -----------------------------------------------------------------------------
# 02_magma_array.sh -- MAGMA gene analysis, one array task per trait.
#
#   sbatch 02_magma_array.sh                 # all 40 traits, 20 at a time
#   sbatch --array=34 02_magma_array.sh      # just SCZ
#
# Array task n maps to DATA row n of traits.tsv, i.e. FILE line n+1 (line 1 is
# the header). SCZ is data row 34 / file line 35. To find a trait's task id:
#   awk 'NR>1 && $1=="SCZ"{print NR-1}' traits.tsv
#
# Per task:
#   1. LDSC sumstats  ->  SNP/P/N table          (02a_sumstats_to_pval.py)
#   2. magma --gene-annot ... --pval ...         -> <trait>.genes.out
#
# IDEMPOTENT: a trait whose .genes.out already exists is skipped, so a partial
# array can simply be resubmitted. Nothing is ever deleted or overwritten.
#
# Outputs:
#   $MAGMA_PVAL/<trait>.pval.tsv   (scratch, ~30 MB each)
#   $MAGMA_OUT/<trait>.genes.out   (project tree, ~1 MB each -- the deliverable)
#
# Runtime is ~5-15 min per trait at 1.0-1.2 M SNPs and ~19 k genes.
#
# RESOURCES: 2 CPU / 8 G. MAGMA's gene analysis is single-threaded -- extra
# cores buy nothing here, so the second core just covers the python Z->P step.
# Peak memory is the ~1.1 M-row P table plus MAGMA's per-gene SNP window, well
# under 8 G. 20 concurrent tasks (rather than 10) because small jobs schedule
# fast on `shared`; 40 traits should clear in ~2 waves.
# -----------------------------------------------------------------------------
set -euo pipefail
# NOTE: under sbatch, SLURM copies this script to a spool directory, so "$0"
# does NOT resolve to this file's real location -- $(dirname $(readlink -f $0))
# yields /var/spool/slurm/... and config.sh is not found there. Hardcode the
# code directory instead. Override with CODE_DIR=... if you relocate the tree.
CODE_DIR=${CODE_DIR:-/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS}
source "$CODE_DIR/config.sh"

echo "**** Job starts ****"; date; hostname
echo "task ${SLURM_ARRAY_TASK_ID:-1}"

module load magma/1.10
activate_env

mkdir -p "$MAGMA_PVAL" "$MAGMA_OUT" "$CODE/logs"

# traits.tsv line 1 is a header, so task k maps to line k+1
LINE=$(( ${SLURM_ARRAY_TASK_ID:-1} + 1 ))
TRAIT=$(awk -v n="$LINE" 'NR==n{print $1}' "$TRAITS_TSV")
SS=$(awk -v n="$LINE" 'NR==n{print $2}' "$TRAITS_TSV")

if [ -z "${TRAIT:-}" ]; then
    echo "no trait on line $LINE of $TRAITS_TSV -- nothing to do"; exit 0
fi
echo "TRAIT=$TRAIT  SUMSTATS=$SS"

# traits.tsv column 2 is normally a bare filename resolved against $GWAS_SRC
# (the shared, read-only gwas_brain directory). An ABSOLUTE path is also
# accepted, which is how traits whose sumstats we downloaded ourselves --
# e.g. PTSD_F3 in $WORKDIR/gwas_new/ -- are pointed at. Nothing is ever
# written back to $GWAS_SRC.
case "$SS" in
    /*) SRC=$SS ;;
    *)  SRC=$GWAS_SRC/$SS ;;
esac
if [ ! -r "$SRC" ]; then
    echo "ERROR: cannot read $SRC"; exit 1
fi

OUT=$MAGMA_OUT/$TRAIT
if [ -f "$OUT.genes.out" ]; then
    echo "$OUT.genes.out already exists -- skipping"
    echo "**** Job ends ****"; date; exit 0
fi

# --- 1. P values -------------------------------------------------------------
PV=$MAGMA_PVAL/$TRAIT.pval.tsv
if [ -f "$PV" ]; then
    echo "[1/2] $PV exists -- reusing"
else
    echo "[1/2] converting Z -> P"
    "$PY" "$CODE/02a_sumstats_to_pval.py" "$SRC" "$PV.tmp$$"
    mv "$PV.tmp$$" "$PV"
fi

NSNP=$(( $(wc -l < "$PV") - 1 ))
echo "  SNPs: $NSNP"

# --- 2. MAGMA gene analysis --------------------------------------------------
# ncol=N takes the per-SNP sample size from the third column of $PV.
echo "[2/2] magma --gene-annot"
magma --bfile "$MAGMA_REF/g1000_eur" \
      --gene-annot "$MAGMA_REF/magma_annot.genes.annot" \
      --pval "$PV" ncol=N \
      --out "$OUT.tmp$$"

for ext in genes.out genes.raw log; do
    [ -f "$OUT.tmp$$.$ext" ] && mv "$OUT.tmp$$.$ext" "$OUT.$ext"
done

echo "---- $TRAIT.genes.out ----"
head -3 "$OUT.genes.out"
echo "genes tested: $(( $(grep -vc '^#' "$OUT.genes.out") - 1 ))"

echo "**** Job ends ****"; date

#!/bin/bash
# -----------------------------------------------------------------------------
# config.sh -- single source of truth for every path in the scDRS pipeline.
# Sourced by all 0*.sh scripts. Nothing here executes work.
#
# SCOPE RULE FOR THIS DIRECTORY:
#   Scripts here only ever CREATE files under
#     - $CODE      (this directory)                 -- code + small text logs
#     - $WORKDIR   (processed-data/.../scDRS)       -- results
#     - $SCRATCH   (/users/mtotty/claude_scratch)   -- heavy intermediates
#   Nothing is ever deleted, moved, or overwritten outside those three.
#   The shared GWAS sumstats and reference panels are read-only inputs and are
#   only ever read or symlinked, never copied or modified.
#
#   >>> NOTE: $WORKDIR and $SCRATCH are OUTSIDE this code directory. They are
#   >>> new, previously non-existent paths that mirror where gsMap already
#   >>> writes its results, and only ever get created -- never overwriting
#   >>> anything. If you would rather keep every output inside the scDRS code
#   >>> folder instead, change WORKDIR/SCRATCH below; nothing else needs edits.
# -----------------------------------------------------------------------------

PROJ=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala
CODE=$PROJ/code/Visium/16_LDSC/scDRS
WORKDIR=$PROJ/processed-data/Visium/16_LDSC/scDRS
SCRATCH=/users/mtotty/claude_scratch/scdrs

# ---- environment -------------------------------------------------------------
CONDA_MODULE=conda/3-24.3.0
ENVDIR=$SCRATCH/envs/scdrs
PY=$ENVDIR/bin/python
SCDRS=$ENVDIR/bin/scdrs

# ---- read-only shared inputs -------------------------------------------------
GWAS_SRC=/dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain
GENE_LOC_SRC=/dcs04/lieber/shared/statsgen/PRS/ref_EUR/NCBI37.3.gene.loc   # hg19 / NCBI37
PLINK_SRC=/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink
PLINK_PREFIX=1000G.EUR.QC          # per-chromosome: ${PLINK_PREFIX}.{1..22}

# spatial inputs -- the SAME h5ad files gsMap consumes
ST_DIR=$PROJ/processed-data/Visium/16_LDSC/gsMap/ST
ANNOT_COL=BS_k16_Semisupervised_wAI

# ---- derived MAGMA paths (all in scratch) -----------------------------------
MAGMA_REF=$SCRATCH/magma_ref               # merged plink panel + snp-loc + annot
MAGMA_PVAL=$SCRATCH/magma_pval             # per-trait SNP/P/N files
MAGMA_OUT=$WORKDIR/magma                   # per-trait .genes.out (small, keep)

# ---- scDRS paths -------------------------------------------------------------
GS_DIR=$WORKDIR/gs                         # munged .gs gene-set files
H5AD_SCDRS=$WORKDIR/h5ad/visium_amygdala_scdrs.h5ad
COV_FILE=$WORKDIR/h5ad/cov.tsv
SCORE_DIR=$WORKDIR/score                   # per-trait .score.gz / .full_score.gz
DOWNSTREAM_DIR=$WORKDIR/downstream         # group-level stats per domain
FIG_DIR=$WORKDIR/figures

# ---- MAGMA gene window (kb up, kb down) -------------------------------------
# 10,1.5 is the MAGMA default and what the scDRS authors use for their gene sets.
MAGMA_WINDOW="10,1.5"

# ---- scDRS parameters --------------------------------------------------------
N_CTRL=1000        # control gene sets
GS_NMAX=1000       # top-N genes retained per trait by munge-gs

# ---- SLURM ------------------------------------------------------------------
PARTITION=shared
WALLTIME=3-00:00:00

TRAITS_TSV=$CODE/traits.tsv

activate_env() {
    module load $CONDA_MODULE
    eval "$(conda shell.bash hook)"
    conda activate $ENVDIR
}

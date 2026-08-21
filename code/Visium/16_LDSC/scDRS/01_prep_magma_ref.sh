#!/bin/bash
#SBATCH --job-name=scdrs_magma_ref
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/01_prep_magma_ref.%j.log
# -----------------------------------------------------------------------------
# 01_prep_magma_ref.sh -- build everything MAGMA needs, ONCE, in scratch.
#
#   sbatch 01_prep_magma_ref.sh
#
# Produces (all under $MAGMA_REF, i.e. scratch -- nothing in the project tree):
#   g1000_eur.{bed,bim,fam}        merged 1000G EUR Phase 3 panel (22 autosomes)
#   snp_loc.txt                    SNP CHR BP, from the merged .bim
#   gene_loc_symbol.txt            NCBI37.3.gene.loc with SYMBOL as the gene ID
#   magma_annot.genes.annot        SNP -> gene assignment, window 10 kb up / 1.5 kb down
#
# WHY THESE INPUTS
#   The plink panel is the same 1000G EUR Phase 3 build the classic s-LDSC
#   pipeline in ../LDSC/ uses, so the LD model is shared across all three
#   methods. NCBI37.3.gene.loc is hg19/NCBI37 -- matching gene_meta_hg19.txt
#   used by ../LDSC/03_bed.R and the gencode.v46lift37 GTF used by gsMap.
#
#   MAGMA writes the col-1 identifier of the gene-loc file into .genes.out.
#   scDRS matches gene sets to the h5ad by SYMBOL, so we put the symbol in
#   col 1 rather than translating Entrez IDs downstream.
#
# Source files are read-only shared data and are only read, never modified.
# -----------------------------------------------------------------------------
set -euo pipefail
# NOTE: under sbatch, SLURM copies this script to a spool directory, so "$0"
# does NOT resolve to this file's real location -- $(dirname $(readlink -f $0))
# yields /var/spool/slurm/... and config.sh is not found there. Hardcode the
# code directory instead. Override with CODE_DIR=... if you relocate the tree.
CODE_DIR=${CODE_DIR:-/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS}
source "$CODE_DIR/config.sh"

echo "**** Job starts ****"; date; hostname

module load magma/1.10
module load plink/1.90b

# GUARD: /usr/bin/plink is PuTTY's SSH client, not genetics plink. The module
# refuses to load on login nodes, so without this check a login-node run would
# silently call the wrong binary. Fail loudly instead.
if ! plink --version 2>&1 | grep -q '^PLINK'; then
    echo "ERROR: 'plink' on PATH is not genetics plink -- got:"
    plink --version 2>&1 | head -3
    echo "Are you on a compute node? plink/1.90b will not load on a login node."
    exit 1
fi
echo "plink: $(plink --version 2>&1 | head -1)"

mkdir -p "$MAGMA_REF" "$CODE/logs"
cd "$MAGMA_REF"

# --- 1. merge the 22 per-chromosome plink files into one panel ---------------
if [ -f "$MAGMA_REF/g1000_eur.bed" ]; then
    echo "[1/4] merged panel already present -- skipping"
else
    echo "[1/4] merging 1000G EUR Phase 3 autosomes"
    : > merge_list.txt
    for chr in $(seq 2 22); do
        echo "$PLINK_SRC/${PLINK_PREFIX}.${chr}" >> merge_list.txt
    done
    plink --bfile "$PLINK_SRC/${PLINK_PREFIX}.1" \
          --merge-list merge_list.txt \
          --make-bed --out g1000_eur
fi
wc -l g1000_eur.bim

# --- 2. SNP location file ----------------------------------------------------
echo "[2/4] writing snp_loc.txt"
awk '{print $2"\t"$1"\t"$4}' g1000_eur.bim > snp_loc.txt
wc -l snp_loc.txt

# --- 3. gene location file keyed on symbol ----------------------------------
# NCBI37.3.gene.loc columns: ENTREZ CHR START STOP STRAND SYMBOL
# Keep autosomes only (LDSC panel is autosomal); drop duplicate symbols,
# keeping the first occurrence so the file stays a valid MAGMA gene-loc.
echo "[3/4] writing gene_loc_symbol.txt"
awk -F'\t' 'BEGIN{OFS="\t"}
     $2 ~ /^([1-9]|1[0-9]|2[0-2])$/ && !seen[$6]++ {print $6,$2,$3,$4,$5}' \
     "$GENE_LOC_SRC" > gene_loc_symbol.txt
wc -l gene_loc_symbol.txt
head -3 gene_loc_symbol.txt

# --- 4. annotate -------------------------------------------------------------
echo "[4/4] magma --annotate  (window $MAGMA_WINDOW)"
magma --annotate window="$MAGMA_WINDOW" \
      --snp-loc snp_loc.txt \
      --gene-loc gene_loc_symbol.txt \
      --out magma_annot

echo "---- annot summary ----"
head -3 magma_annot.genes.annot
echo "genes annotated: $(( $(wc -l < magma_annot.genes.annot) - 2 ))"

echo "**** Job ends ****"; date

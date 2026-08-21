#!/bin/bash
# =============================================================================
# 00_download_resources.sh -- fetch + unpack the gsMap reference resource bundle
#
# Project : spatialAmygdala (LIBD) -- Visium / 16_LDSC / gsMap
# Purpose : download gsMap_resource.tar.gz (tens of GB) from the Yang Lab, extract
#           it into scratch, and write a manifest of the exact absolute paths the
#           downstream gsMap scripts need.
#
# Where it runs : COMPUTE NODE ONLY (never the login node -- this is a large
#                 network transfer plus a multi-GB tar extraction).
#                   srun --pty -p shared --cpus-per-task 4 --mem 16G --time 1-00:00:00 bash
#
# Idempotent : re-running skips the download if the tarball is already complete
#              (wget -c resumes), and skips extraction if the sentinel file
#              quick_mode/snp_gene_weight_matrix.h5ad already exists.
#              Nothing is ever deleted by this script.
# =============================================================================
#SBATCH --job-name=gsmap_resources
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00

set -euo pipefail

# ---- hardcoded paths --------------------------------------------------------
SCRATCH=/users/mtotty/claude_scratch/gsmap
URL=https://yanglab.westlake.edu.cn/data/gsMap/gsMap_resource.tar.gz
TARBALL=${SCRATCH}/gsMap_resource.tar.gz
RES_DIR=${SCRATCH}/gsMap_resource
SENTINEL=${RES_DIR}/quick_mode/snp_gene_weight_matrix.h5ad
MANIFEST=${SCRATCH}/gsMap_resource_manifest.txt

mkdir -p "${SCRATCH}"
cd "${SCRATCH}"

# ---- 1. download ------------------------------------------------------------
if [[ -f "${SENTINEL}" ]]; then
  echo "[resources] sentinel present (${SENTINEL}) -- skipping download + extract"
else
  echo "[resources] downloading ${URL}"
  # -c resumes a partial transfer; --tries/--timeout keep a stalled socket from
  # burning the whole walltime.
  wget -c --tries=10 --timeout=60 --progress=dot:giga -O "${TARBALL}" "${URL}"
  ls -lh "${TARBALL}"

  # ---- 2. extract -----------------------------------------------------------
  echo "[resources] extracting into ${SCRATCH}"
  tar -xzvf "${TARBALL}" -C "${SCRATCH}" | tail -50
fi

# ---- 3. verify the expected tree -------------------------------------------
echo "[resources] verifying tree"
MISSING=0
for P in \
  "${RES_DIR}/genome_annotation/gtf" \
  "${RES_DIR}/genome_annotation/enhancer" \
  "${RES_DIR}/LD_Reference_Panel/1000G_EUR_Phase3_plink" \
  "${RES_DIR}/LDSC_resource/hapmap3_snps" \
  "${RES_DIR}/LDSC_resource/weights_hm3_no_hla" \
  "${RES_DIR}/quick_mode/baseline" \
  "${RES_DIR}/quick_mode/SNP_gene_pair" \
  "${RES_DIR}/quick_mode/snp_gene_weight_matrix.h5ad" ; do
  if [[ -e "${P}" ]]; then echo "  OK      ${P}"; else echo "  MISSING ${P}"; MISSING=1; fi
done

# ---- 4. manifest ------------------------------------------------------------
echo "[resources] writing ${MANIFEST}"
{
  echo "===== gsMap resource bundle manifest ====="
  echo "date        : $(date)"
  echo "host        : $(hostname)"
  echo "source url  : ${URL}"
  echo "tarball     : ${TARBALL}"
  echo "resource dir: ${RES_DIR}"
  echo

  echo "----- top-level tree (2 levels) -----"
  find "${RES_DIR}" -maxdepth 2 -type d | sort
  echo

  echo "----- du -sh of each top-level subdir -----"
  for D in "${RES_DIR}"/*/ ; do du -sh "${D}"; done
  echo
  echo "----- du -sh total -----"
  du -sh "${RES_DIR}"
  echo

  echo "----- genome_annotation/gtf contents (GTF filename) -----"
  ls -lh "${RES_DIR}/genome_annotation/gtf/" 2>&1
  echo
  echo "----- genome_annotation/enhancer contents -----"
  ls -lh "${RES_DIR}/genome_annotation/enhancer/" 2>&1
  echo

  echo "----- LD_Reference_Panel/1000G_EUR_Phase3_plink (first 12 files; bfile prefix) -----"
  ls -lh "${RES_DIR}/LD_Reference_Panel/1000G_EUR_Phase3_plink/" 2>&1 | head -14
  echo "bed/bim/fam count: $(ls "${RES_DIR}/LD_Reference_Panel/1000G_EUR_Phase3_plink/" | grep -cE '\.(bed|bim|fam)$' || true)"
  echo

  echo "----- LDSC_resource/hapmap3_snps (first 8) -----"
  ls -lh "${RES_DIR}/LDSC_resource/hapmap3_snps/" 2>&1 | head -10
  echo
  echo "----- LDSC_resource/weights_hm3_no_hla (first 8) -----"
  ls -lh "${RES_DIR}/LDSC_resource/weights_hm3_no_hla/" 2>&1 | head -10
  echo

  echo "----- quick_mode/ contents -----"
  ls -lh "${RES_DIR}/quick_mode/" 2>&1
  echo
  echo "----- quick_mode/baseline (first 10) -----"
  ls -lh "${RES_DIR}/quick_mode/baseline/" 2>&1 | head -12
  echo
  echo "----- quick_mode/SNP_gene_pair (first 10) -----"
  ls -lh "${RES_DIR}/quick_mode/SNP_gene_pair/" 2>&1 | head -12
  echo

  echo "----- homologs/ -----"
  ls -lh "${RES_DIR}/homologs/" 2>&1 | head -10
  echo

  echo "----- GTF genome-build check (head of GTF) -----"
  GTF=$(ls "${RES_DIR}/genome_annotation/gtf/"*.gtf* 2>/dev/null | head -1 || true)
  echo "GTF file: ${GTF}"
  if [[ -n "${GTF}" ]]; then
    case "${GTF}" in
      *.gz) zcat "${GTF}" | head -8 ;;
      *)    head -8 "${GTF}" ;;
    esac
  fi
  echo

  echo "===== PATHS FOR DOWNSTREAM gsMap SCRIPTS ====="
  echo "--gsMap_resource_dir  ${RES_DIR}"
  echo "--bfile_root          ${RES_DIR}/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC"
  echo "--keep_snp_root       ${RES_DIR}/LDSC_resource/hapmap3_snps/hm"
  echo "--gtf_annotation_file ${GTF}"
  echo "--w_file              ${RES_DIR}/LDSC_resource/weights_hm3_no_hla/weights."
  echo "(verify the trailing prefix stems against the ls listings above)"
  echo
  echo "missing_components_flag: ${MISSING}"
} > "${MANIFEST}"

echo "[resources] DONE (missing flag = ${MISSING})"
echo "[resources] manifest: ${MANIFEST}"

# =============================================================================
# 01_spe_to_h5ad.R
#
# Export the BayesSpace k16 semisupervised SPE to one h5ad per capture area
# (stitched section = colData$sample_id) for gsMap.
#
# gsMap input contract satisfied here:
#   * raw integer counts in  adata.layers['count']   (gsMap --data_layer count)
#     (the same matrix is also left in adata.X for convenience)
#   * spatial coordinates in adata.obsm['spatial']   (numeric, 2 columns)
#   * annotation in          adata.obs['BS_k16_Semisupervised_wAI']
#                                                    (gsMap --annotation)
#   * var_names = gene SYMBOLS -- gsMap matches genes by symbol
#   * species is HUMAN, so no --homolog_file is needed downstream
#
# Gene universe follows the classic stratified-LDSC pipeline convention
# (code/Visium/16_LDSC/LDSC/01_aggregate.R): protein_coding only, then
# deduplicated on rowData$gene_name.
#
# NOTE ON SPLIT UNIT: this SPE is a visiumStitched object. colData$sample_id
# has 7 levels (donors) and colData$capture_area has 51. spatialCoords() holds
# the STITCHED fullres coordinates, which are only mutually consistent within a
# sample_id; the per-capture-area coordinates are kept separately as
# pxl_{row,col}_in_fullres_original. We therefore split on sample_id, so each
# h5ad is one contiguous spatial frame -- which is what gsMap's spatial KNN
# graph needs.
#
# Run under `module load conda_R/4.5` on a COMPUTE NODE (login nodes kill R).
# This script writes ONLY into
#   processed-data/Visium/16_LDSC/gsMap/ST/   and   code/Visium/16_LDSC/gsMap/
# It deletes nothing.
# =============================================================================

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(zellkonverter)
})

## ---- hardcoded paths (repo deliverable) -----------------------------------
proj_dir   <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
spe_rds    <- file.path(proj_dir, "processed-data", "Visium", "07_clustering",
                        "BayesSpace", "MarkerGenes",
                        "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds")
gsmap_code <- file.path(proj_dir, "code", "Visium", "16_LDSC", "gsMap")
workdir    <- file.path(proj_dir, "processed-data", "Visium", "16_LDSC", "gsMap")
st_dir     <- file.path(workdir, "ST")

manifest   <- file.path(gsmap_code, "01_h5ad_manifest.txt")

## annotation column = the gsMap --annotation
ann_col     <- "BS_k16_Semisupervised_wAI"
## split unit = one stitched section per donor
sample_col  <- "sample_id"
## 01_aggregate.R convention: drop domains with fewer than this many spots
min_spots_per_domain <- 70

dir.create(st_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
cat("reading:", spe_rds, "\n")
spe <- readRDS(spe_rds)
cat("loaded SPE:", nrow(spe), "genes x", ncol(spe), "spots\n")

## ---- 1. drop domains with <70 spots GLOBALLY -------------------------------
## (mirrors code/Visium/16_LDSC/LDSC/01_aggregate.R)
dom_counts   <- table(colData(spe)[[ann_col]])
dom_to_drop  <- names(dom_counts[dom_counts < min_spots_per_domain])
cat("domains dropped (<", min_spots_per_domain, " spots globally): ",
    if (length(dom_to_drop)) paste(dom_to_drop, collapse = ", ") else "none",
    "\n", sep = "")
spe <- spe[, !(colData(spe)[[ann_col]] %in% dom_to_drop)]
colData(spe)[[ann_col]] <- factor(as.character(colData(spe)[[ann_col]]))
domains_kept <- levels(colData(spe)[[ann_col]])
cat("domains retained (", length(domains_kept), "): ",
    paste(domains_kept, collapse = ", "), "\n", sep = "")

## ---- 2. gene universe: protein_coding, deduplicated on gene_name -----------
stopifnot(all(c("gene_name", "gene_type") %in% names(rowData(spe))))
keep_pc  <- rowData(spe)$gene_type == "protein_coding"
spe      <- spe[keep_pc, ]
cat("protein_coding genes:", nrow(spe), "\n")
spe      <- spe[!duplicated(rowData(spe)$gene_name), ]
cat("after dedup on gene_name:", nrow(spe), "\n")

## var_names must be the SYMBOL -- gsMap matches genes by symbol
rownames(spe) <- as.character(rowData(spe)$gene_name)
stopifnot(!any(is.na(rownames(spe))), !any(duplicated(rownames(spe))))
n_genes <- nrow(spe)

## ---- 3. build a tidy per-spot annotation table -----------------------------
cd  <- colData(spe)
obs <- DataFrame(row.names = colnames(spe))
obs[[ann_col]]  <- factor(as.character(cd[[ann_col]]), levels = domains_kept)
obs$sample_id   <- factor(as.character(cd[[sample_col]]))
## carried through for traceability / optional downstream filtering:
if ("capture_area" %in% names(cd))
  obs$capture_area <- factor(as.character(cd$capture_area))
if ("exclude_overlapping" %in% names(cd))
  obs$exclude_overlapping <- as.logical(cd$exclude_overlapping)
if ("in_tissue" %in% names(cd)) obs$in_tissue <- as.logical(cd$in_tissue)
if ("sum_umi"  %in% names(cd))  obs$sum_umi   <- as.numeric(cd$sum_umi)
if ("sum_gene" %in% names(cd))  obs$sum_gene  <- as.numeric(cd$sum_gene)
if ("key"      %in% names(cd))  obs$key       <- as.character(cd$key)

## hard requirement: gsMap cannot have NA annotation
stopifnot(!any(is.na(obs[[ann_col]])))

## ---- 4. spatial coordinates -------------------------------------------------
coords <- spatialCoords(spe)
storage.mode(coords) <- "double"
stopifnot(ncol(coords) == 2, !any(is.na(coords)))
cat("spatialCoords columns:", paste(colnames(coords), collapse = ", "), "\n")

## ---- 5. counts --------------------------------------------------------------
cnt <- assay(spe, "counts")
cat("counts class:", class(cnt)[1], "\n")

## sanity check: counts must be raw integers for gsMap
nz <- cnt@x
cat("counts nonzero range:", min(nz), "-", max(nz),
    "| all integer-valued:", all(nz == floor(nz)), "\n")
stopifnot(all(nz == floor(nz)), min(nz) >= 0)
rm(nz); invisible(gc())

## ---- 6. optional: record the classic-LDSC gene universe for comparison ------
## 01_aggregate.R additionally applies edgeR::filterByExpr on the domain
## pseudobulk, landing on 16440 genes. We do NOT apply it to the h5ad (gsMap
## does its own gene ranking), but we write the list so the two pipelines'
## gene universes can be reconciled later.
gene_universe_file <- file.path(gsmap_code, "01_gene_universe_classic.txt")
invisible(try({
  suppressPackageStartupMessages({ library(scuttle); library(edgeR) })
  sce_tmp <- SingleCellExperiment(list(counts = cnt), colData = obs)
  pb <- scuttle::aggregateAcrossCells(sce_tmp, ids = obs[[ann_col]])
  hi <- edgeR::filterByExpr(pb, group = pb[[ann_col]])
  writeLines(sort(rownames(pb)[hi]), gene_universe_file)
  cat("classic filterByExpr gene universe:", sum(hi), "genes ->",
      gene_universe_file, "\n")
  rm(sce_tmp, pb); invisible(gc())
}, silent = FALSE))

## ---- 7. split by sample_id and write one h5ad each --------------------------
samples <- sort(unique(as.character(obs$sample_id)))
cat("\ncapture areas (stitched sections) to export:",
    paste(samples, collapse = ", "), "\n\n")

mf <- file(manifest, open = "wt")
writeLines(c(
  "=============================================================",
  "gsMap spatial input manifest -- per-capture-area h5ad",
  paste("generated:", format(Sys.time())),
  paste("source SPE:", spe_rds),
  paste("output dir:", st_dir),
  "",
  paste("split column      :", sample_col, "(visiumStitched section)"),
  paste("annotation column :", ann_col),
  paste("gene universe     : protein_coding, deduplicated on gene_name"),
  paste("n genes           :", n_genes),
  paste("domain filter     : >=", min_spots_per_domain, "spots globally"),
  paste("domains dropped   :",
        if (length(dom_to_drop)) paste(dom_to_drop, collapse = ", ") else "none"),
  paste("domains retained  :", length(domains_kept), "-",
        paste(domains_kept, collapse = ", ")),
  "",
  "layout of each h5ad:",
  "  X                 = raw counts (same matrix as layers['count'])",
  "  layers['count']   = raw integer counts   (gsMap --data_layer count)",
  "  obsm['spatial']   = stitched fullres pixel coords (n_spots x 2)",
  paste0("  obs['", ann_col, "'] = domain label (no NA)"),
  "  var_names         = gene symbols",
  "============================================================="
), mf)

written <- character(0)
for (s in samples) {
  idx <- which(as.character(obs$sample_id) == s)
  cat("---", s, ":", length(idx), "spots\n")

  m_s   <- cnt[, idx, drop = FALSE]
  obs_s <- obs[idx, , drop = FALSE]
  crd_s <- coords[idx, , drop = FALSE]

  ## obs_names must be unique within the file; barcodes repeat across the
  ## capture areas that make up a stitched section, so prefer `key`.
  nm <- if (!is.null(obs_s$key) && !any(duplicated(obs_s$key))) {
    as.character(obs_s$key)
  } else {
    make.unique(paste(s, colnames(m_s), sep = "_"))
  }
  colnames(m_s) <- nm
  rownames(obs_s) <- nm
  rownames(crd_s) <- nm

  ## drop factor levels that this section does not contain, but keep the
  ## annotation as a factor so it lands as a categorical in obs
  obs_s$sample_id <- factor(as.character(obs_s$sample_id))
  if (!is.null(obs_s$capture_area))
    obs_s$capture_area <- factor(as.character(obs_s$capture_area))

  ## two assays pointing at the same matrix: X_name="counts" sends one to X,
  ## the other becomes layers[['count']] -- the layer gsMap reads.
  sce_s <- SingleCellExperiment(
    assays  = list(counts = m_s, count = m_s),
    colData = obs_s,
    rowData = DataFrame(
      gene_name = rownames(m_s),
      gene_id   = as.character(rowData(spe)$gene_id),
      gene_type = as.character(rowData(spe)$gene_type),
      row.names = rownames(m_s)
    )
  )
  ## reducedDims -> obsm, so this becomes adata.obsm['spatial'].
  ## Strip dimnames first: if the matrix carries row/col names, zellkonverter
  ## writes obsm['spatial'] as a structured array and anndata reads it back as
  ## a pandas DataFrame, which breaks gsMap's positional coords[:, 0] indexing.
  ## Without dimnames it round-trips as a plain float64 ndarray.
  crd_plain <- matrix(as.numeric(crd_s), ncol = 2, dimnames = NULL)
  reducedDim(sce_s, "spatial") <- crd_plain

  out_h5 <- file.path(st_dir, paste0(s, ".h5ad"))
  writeH5AD(sce_s, file = out_h5, X_name = "counts", compression = "gzip")
  written <- c(written, out_h5)

  ## ---- manifest entry ----
  dt <- table(as.character(obs_s[[ann_col]]))
  dt <- dt[order(-dt)]
  writeLines(c(
    "",
    paste0("### ", s),
    paste("file       :", out_h5),
    paste("size       :", format(structure(file.size(out_h5),
                                           class = "object_size"),
                                 units = "auto")),
    paste("n_spots    :", ncol(m_s)),
    paste("n_genes    :", nrow(m_s)),
    paste("n_capture_areas:",
          if (!is.null(obs_s$capture_area)) nlevels(obs_s$capture_area) else NA),
    paste("counts sum :", format(sum(m_s), scientific = FALSE)),
    paste("n exclude_overlapping TRUE:",
          if (!is.null(obs_s$exclude_overlapping))
            sum(obs_s$exclude_overlapping, na.rm = TRUE) else NA),
    paste("coord range: x [", paste(range(crd_s[, 1]), collapse = ", "),
          "]  y [", paste(range(crd_s[, 2]), collapse = ", "), "]"),
    paste("domains present:", length(dt), "of", length(domains_kept)),
    paste("  missing      :",
          paste(setdiff(domains_kept, names(dt)), collapse = ", ")),
    "  domain spot counts:",
    paste0("    ", names(dt), ": ", as.integer(dt))
  ), mf)
  flush(mf)

  rm(m_s, obs_s, crd_s, crd_plain, sce_s); invisible(gc())
}

writeLines(c("", "=============================================================",
             paste("wrote", length(written), "h5ad files"),
             written, "=== END OF MANIFEST ==="), mf)
close(mf)

cat("\nwrote", length(written), "h5ad files into", st_dir, "\n")
print(written)
cat("manifest:", manifest, "\n")
cat("DONE\n")

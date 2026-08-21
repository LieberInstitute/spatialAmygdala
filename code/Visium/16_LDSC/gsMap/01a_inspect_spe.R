# =============================================================================
# 01a_inspect_spe.R
# Read-only inspection of the BayesSpace k16 semisupervised SPE, to plan the
# gsMap spatial-input export (per-capture-area h5ad).
#
# Writes a plain-text report to 01a_spe_summary.txt in this directory.
# Run under `module load conda_R/4.5` on a COMPUTE NODE (login nodes kill R).
#
# This script only READS the SPE. It creates no files outside 16_LDSC/gsMap/.
# =============================================================================

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SingleCellExperiment)
})

## ---- hardcoded paths (repo deliverable) -----------------------------------
proj_dir   <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
spe_rds    <- file.path(proj_dir, "processed-data", "Visium", "07_clustering",
                        "BayesSpace", "MarkerGenes",
                        "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds")
gsmap_code <- file.path(proj_dir, "code", "Visium", "16_LDSC", "gsMap")
out_txt    <- file.path(gsmap_code, "01a_spe_summary.txt")

## annotation column used as the gsMap --annotation
ann_col <- "BS_k16_Semisupervised_wAI"

cat("reading:", spe_rds, "\n")
spe <- readRDS(spe_rds)

## ---- everything below is captured into the report -------------------------
con <- file(out_txt, open = "wt")
sink(con, split = TRUE)

cat("=============================================================\n")
cat("SPE inspection report for gsMap input prep\n")
cat("generated:", format(Sys.time()), "\n")
cat("source rds:", spe_rds, "\n")
cat("=============================================================\n\n")

cat("### class\n"); print(class(spe)); cat("\n")

cat("### dim(spe)  [genes x spots]\n"); print(dim(spe)); cat("\n")

cat("### assayNames(spe)\n"); print(assayNames(spe)); cat("\n")

cat("### names(colData(spe))\n"); print(names(colData(spe))); cat("\n")

## ---- locate the sample id column ------------------------------------------
cat("### candidate sample-id columns\n")
cand <- grep("sample|slide|capture|brnum|section|array",
             names(colData(spe)), ignore.case = TRUE, value = TRUE)
print(cand)
for (cc in cand) {
  nu <- length(unique(colData(spe)[[cc]]))
  cat(sprintf("  %-30s n_unique = %d\n", cc, nu))
}
cat("\n")

## prefer literal "sample_id" (SpatialExperiment standard); else first candidate
sample_col <- if ("sample_id" %in% names(colData(spe))) "sample_id" else cand[1]
cat("### chosen sample column:", sample_col, "\n\n")

sid <- as.character(colData(spe)[[sample_col]])

cat("### spots per capture area (", sample_col, ")\n", sep = "")
tab_s <- sort(table(sid), decreasing = TRUE)
print(tab_s)
cat("n capture areas =", length(tab_s), "\n")
cat("total spots     =", sum(tab_s), "\n\n")

## ---- annotation column -----------------------------------------------------
cat("### annotation column:", ann_col, "\n")
cat("class:\n"); print(class(colData(spe)[[ann_col]]))
ann <- as.character(colData(spe)[[ann_col]])
cat("n NA in annotation =", sum(is.na(ann)), "\n\n")

cat("### global spots per domain (", ann_col, "), sorted\n", sep = "")
tab_d <- sort(table(ann), decreasing = TRUE)
print(tab_d)
cat("n domains =", length(tab_d), "\n\n")

cat("### domains with <70 spots globally (dropped by 01_aggregate.R convention)\n")
print(names(tab_d)[tab_d < 70])
cat("\n")
cat("### domains RETAINED (>=70 spots globally)\n")
print(sort(names(tab_d)[tab_d >= 70]))
cat("n retained =", sum(tab_d >= 70), "\n\n")

cat("### table(sample_id, ", ann_col, ")\n", sep = "")
xt <- table(sid, ann)
print(xt)
cat("\n")

cat("### same table restricted to RETAINED domains (>=70 globally)\n")
keep_dom <- names(tab_d)[tab_d >= 70]
xt_keep <- xt[, colnames(xt) %in% keep_dom, drop = FALSE]
print(xt_keep)
cat("\n")

cat("### per-sample coverage of retained domains\n")
cat("(n_domains_present = retained domains with >0 spots;\n")
cat(" n_domains_ge10    = retained domains with >=10 spots)\n")
cov <- data.frame(
  sample_id         = rownames(xt_keep),
  n_spots           = as.integer(rowSums(xt_keep)),
  n_domains_present = as.integer(rowSums(xt_keep > 0)),
  n_domains_ge10    = as.integer(rowSums(xt_keep >= 10)),
  min_domain_n      = as.integer(apply(xt_keep, 1, min)),
  stringsAsFactors  = FALSE
)
cov <- cov[order(-cov$n_spots), ]
print(cov, row.names = FALSE)
cat("\n")

cat("### missing retained domains per sample\n")
for (s in rownames(xt_keep)) {
  miss <- colnames(xt_keep)[xt_keep[s, ] == 0]
  cat(sprintf("  %-20s missing(%d): %s\n", s, length(miss),
              if (length(miss)) paste(miss, collapse = ", ") else "-"))
}
cat("\n")

## ---- spatial coordinates ----------------------------------------------------
cat("### colnames(spatialCoords(spe))\n")
sc <- spatialCoords(spe)
print(colnames(sc))
cat("class:\n"); print(class(sc))
cat("storage.mode:", storage.mode(sc), "\n")
cat("dim:\n"); print(dim(sc))
cat("head:\n"); print(head(sc, 3))
cat("n NA in coords =", sum(is.na(sc)), "\n\n")

## ---- rowData ----------------------------------------------------------------
cat("### names(rowData(spe))\n"); print(names(rowData(spe))); cat("\n")
cat("### head(rowData(spe), 5)\n")
print(head(as.data.frame(rowData(spe)), 5))
cat("\n")

if ("gene_type" %in% names(rowData(spe))) {
  cat("### table(rowData$gene_type) top 15\n")
  print(head(sort(table(rowData(spe)$gene_type), decreasing = TRUE), 15))
  cat("n protein_coding =", sum(rowData(spe)$gene_type == "protein_coding"), "\n\n")
} else cat("### NO gene_type column in rowData!\n\n")

if ("gene_name" %in% names(rowData(spe))) {
  gn <- rowData(spe)$gene_name
  cat("### gene_name summary\n")
  cat("n genes            =", length(gn), "\n")
  cat("n unique gene_name =", length(unique(gn)), "\n")
  cat("n NA gene_name     =", sum(is.na(gn)), "\n")
  cat("n duplicated       =", sum(duplicated(gn)), "\n")
  if ("gene_type" %in% names(rowData(spe))) {
    pc <- rowData(spe)$gene_type == "protein_coding"
    cat("protein_coding genes                  =", sum(pc), "\n")
    cat("protein_coding after dedup(gene_name) =",
        sum(!duplicated(gn[pc])), "\n")
  }
  cat("head:\n"); print(head(gn, 8)); cat("\n")
} else cat("### NO gene_name column in rowData!\n\n")

cat("### head(rownames(spe))\n"); print(head(rownames(spe), 5)); cat("\n")

## ---- counts assay ------------------------------------------------------------
cat("### counts assay\n")
if ("counts" %in% assayNames(spe)) {
  cnt <- assay(spe, "counts")
  cat("class:\n"); print(class(cnt))
  cat("dim:\n"); print(dim(cnt))
  sub <- cnt[seq_len(min(2000, nrow(cnt))), seq_len(min(2000, ncol(cnt)))]
  sv <- as.numeric(sub@x)
  cat("range of nonzero values in a 2000x2000 corner:\n")
  print(range(sv))
  cat("all integer-valued in that corner:", all(sv == floor(sv)), "\n")
  cat("total sum of that corner:", sum(sv), "\n")
} else cat("NO 'counts' assay -- assays are:", paste(assayNames(spe), collapse=", "), "\n")
cat("\n")

## ---- zellkonverter smoke test -------------------------------------------------
## gsMap needs an h5ad; confirm we can actually write one under this R module
## before committing to the full export in 01_spe_to_h5ad.R.
cat("### zellkonverter availability / smoke test\n")
ok <- requireNamespace("zellkonverter", quietly = TRUE)
cat("zellkonverter installed:", ok, "\n")
if (ok) {
  cat("version:", as.character(packageVersion("zellkonverter")), "\n")
  tmp_h5 <- file.path(gsmap_code, "01a_zellkonverter_smoketest.h5ad")
  res <- try({
    m <- matrix(as.integer(rpois(200, 3)), nrow = 20, ncol = 10)
    rownames(m) <- paste0("G", 1:20); colnames(m) <- paste0("C", 1:10)
    sce_t <- SingleCellExperiment::SingleCellExperiment(list(counts = m))
    SummarizedExperiment::colData(sce_t)$grp <- rep(c("a","b"), 5)
    SingleCellExperiment::reducedDim(sce_t, "spatial") <-
      matrix(as.numeric(1:20), ncol = 2)
    zellkonverter::writeH5AD(sce_t, tmp_h5, X_name = "counts")
    file.size(tmp_h5)
  }, silent = TRUE)
  if (inherits(res, "try-error")) {
    cat("SMOKE TEST FAILED:\n"); cat(as.character(res), "\n")
  } else {
    cat("SMOKE TEST OK -- wrote", res, "bytes to", tmp_h5, "\n")
  }
}
cat("\n")

cat("### sessionInfo()\n")
print(sessionInfo())

cat("\n=== END OF REPORT ===\n")
sink()
close(con)
cat("wrote report to", out_txt, "\n")

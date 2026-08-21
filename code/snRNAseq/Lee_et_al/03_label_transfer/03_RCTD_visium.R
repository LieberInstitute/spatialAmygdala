library(here)
library(SpatialExperiment)
library(spacexr)
library(SingleCellExperiment)
library(Matrix)

plots_dir     <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

LABEL_COL <- "subtype"

# --------------------------------------------------------------------------
# Spatial: SPE -> SpatialRNA
# --------------------------------------------------------------------------
spe <- readRDS(here("processed-data","Visium","07_clustering","BayesSpace",
                    "MarkerGenes","spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
spe <- spe[!duplicated(rownames(spe)), ]

spatial_counts <- counts(spe)
stopifnot(!any(duplicated(rownames(spatial_counts))))

coords <- as.data.frame(spatialCoords(spe))
stopifnot(ncol(coords) == 2)
colnames(coords) <- c("x", "y")
barcodes <- colnames(spatial_counts)
rownames(coords) <- barcodes

nUMI_spatial <- Matrix::colSums(spatial_counts)
names(nUMI_spatial) <- barcodes
stopifnot(identical(rownames(coords), colnames(spatial_counts)),
          identical(names(nUMI_spatial), colnames(spatial_counts)))

puck <- SpatialRNA(coords, spatial_counts, nUMI_spatial)

# --------------------------------------------------------------------------
# Reference: Girgenti SCE -> Reference
# Girgenti is lognorm-only; recover approximate integer counts via expm1.
# --------------------------------------------------------------------------
sce <- readRDS(here("processed-data","snRNAseq","Lee_at_al","girgenti_sce.rds"))

# recover counts from lognorm (same approach as the Yu script)
if (!"counts" %in% assayNames(sce)) {
  # use logcounts (or whichever lognorm assay exists)
  lognorm_assay <- intersect(c("logcounts","X","data"), assayNames(sce))[1]
  counts(sce) <- round(expm1(assay(sce, lognorm_assay)))
} else {
  cmax <- max(counts(sce)[1:min(1000,nrow(sce)), 1:min(50,ncol(sce))])
  if (!isTRUE(all(counts(sce)[1:1000,1:50] == round(counts(sce)[1:1000,1:50]))))
    counts(sce) <- round(expm1(counts(sce)))   # counts slot is actually lognorm
}

sce <- sce[!duplicated(rownames(sce)), ]

# label column + drop tiny cell types (<25 cells)
stopifnot(LABEL_COL %in% colnames(colData(sce)))
sce <- sce[, !is.na(colData(sce)[[LABEL_COL]])]
ct_counts  <- table(colData(sce)[[LABEL_COL]])
keep_types <- names(ct_counts)[ct_counts >= 25]
sce <- sce[, colData(sce)[[LABEL_COL]] %in% keep_types]

ref_counts <- counts(sce)
cell_types <- factor(colData(sce)[[LABEL_COL]])
names(cell_types) <- colnames(sce)
nUMI_ref <- Matrix::colSums(ref_counts)
names(nUMI_ref) <- colnames(sce)

reference <- Reference(ref_counts, cell_types, nUMI_ref)

# ========== RCTD ==========
myRCTD <- create.RCTD(puck, reference, MAX_MULTI_TYPES = 5, max_cores = 30)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "multi")

saveRDS(myRCTD, here(processed_dir, paste0("rctd_girgenti_subtype_", LABEL_COL, ".rds")))
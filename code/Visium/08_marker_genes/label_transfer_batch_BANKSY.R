args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Please provide exactly one sample_id as argument")
query_id <- args[1]

library(Seurat)
library(SpatialExperiment)
library(Matrix)
library(scuttle)
library(here)
library(escheR)
library(SingleCellExperiment)

ref_id <- "Br8325"
label_col <- "BS_manual_ITC_smoothed"

message("Running label transfer: ", ref_id, " → ", query_id)

options(future.globals.maxSize = 2000 * 1024^2)

# ─── Load data ─────────────────────────────────────────
load(here("processed-data", "Visium", "06_batch_correction", "spe_harmony.Rdata"))  # loads `spe`
spe.8325 <- readRDS(here("processed-data", "Visium", "08_marker_genes", "spe_bs_8325_ITC.rds"))
svg_file <- read.csv(here("processed-data", "Visium", "04_feature_selection", "nnSVG_summary.csv"))
top_svg <- svg_file$gene_name[1:2000]

if (!is.null(rowData(spe)$gene_name)) {
    rownames(spe) <- rowData(spe)$gene_name
}

spe.ref <- spe.8325
spe.query <- spe[, spe$sample_id == "Br9280"]


# ─── Subset and normalize ──────────────────────────────
shared_genes <- intersect(top_svg, rownames(spe.ref))
stopifnot(length(shared_genes) >= 1000)


spe.ref <- spe.ref[shared_genes, ]
spe.query <- spe.query[shared_genes, ]


# ─── Combine and Run BANKSY ────────────────────────────
message("Running BANKSY embedding...")


# Run BANKSY
set.seed(1000)
lambda <- 0.4
k_geom <- 18
npcs <- 20
aname <- "logcounts"


spe.ref <- Banksy::computeBanksy(spe.ref, assay_name = aname, k_geom = k_geom)
spe.ref <- Banksy::runBanksyPCA(spe.ref, lambda = lambda, npcs = npcs, group = "capture_area")

spe.query <- Banksy::computeBanksy(spe.query, assay_name = aname, k_geom = k_geom)
spe.query <- Banksy::runBanksyPCA(spe.query, lambda = lambda, npcs = npcs, group = "capture_area")


# ─── Label transfer ────────────────────────────────────
message("Finding anchors and transferring labels...")


seu.ref <- as.Seurat(spe.ref, counts = "counts")
seu.qry <- as.Seurat(spe.query, counts = "counts")


VariableFeatures(seu.ref) <- shared_genes
VariableFeatures(seu.qry) <- shared_genes


seu.ref[[label_col]] <- colData(spe.ref)[[label_col]]
ref_labels <- as.character(seu.ref@meta.data[[label_col]])


anchors <- FindTransferAnchors(
reference = seu.ref,
query = seu.qry,
reduction = "pcaproject",
features = shared_genes,
reference.reduction = "PCA_M0_lam0.4",
dims = 1:npcs,
normalization.method = "SCT"
)


pred <- TransferData(
anchorset = anchors,
refdata = ref_labels,
dims = 1:npcs,
weight.reduction = seu.qry[["PCA_M0_lam0.4"]]
)


seu.qry <- AddMetaData(seu.qry, metadata = pred)
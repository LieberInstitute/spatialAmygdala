args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Please provide exactly one sample_id as argument")
query_id <- args[1]

library(Seurat)
library(SpatialExperiment)
library(Matrix)
library(scuttle)
library(here)
library(escheR)
library(PRECAST)
library(purrr)

ref_id <- "Br8325"
label_col <- "BS_manual_ITC_smoothed"

message("Running label transfer: ", ref_id, " → ", query_id)

# ─── Load data ─────────────────────────────────────────
load(here("processed-data", "Visium", "06_batch_correction", "spe_harmony.Rdata"))
spe.8325 <- readRDS(here("processed-data", "Visium", "08_marker_genes", "spe_bs_8325_ITC.rds"))
svg_file <- read.csv(here("processed-data", "Visium", "04_feature_selection", "nnSVG_summary.csv"))
top_svg <- svg_file$gene_name[1:2000]

if (!is.null(rowData(spe)$gene_name)) rownames(spe) <- rowData(spe)$gene_name

spe.ref <- spe.8325
spe.query <- spe[, spe$sample_id == query_id]

# ─── Subset ────────────────────────────────────────────
shared_genes <- intersect(top_svg, rownames(spe.ref))
stopifnot(length(shared_genes) >= 1000)

# Add spatial coords for PRECAST
colData(spe.ref)$row <- spe.ref$array_row
colData(spe.ref)$col <- spe.ref$array_col
colData(spe.query)$row <- spe.query$array_row
colData(spe.query)$col <- spe.query$array_col

# ─── Build Seurat objects ──────────────────────────────
make_seurat <- function(spe_obj) {
  CreateSeuratObject(
    counts = as.matrix(assay(spe_obj[shared_genes, ], "counts")),
    meta.data = as.data.frame(colData(spe_obj)),
    project = spe_obj$sample_id[1]
  )
}

seu.ref <- make_seurat(spe.ref)
seu.qry <- make_seurat(spe.query)

# ─── Run PRECAST separately ────────────────────────────
run_precast <- function(seu) {
  preobj <- CreatePRECASTObject(seuList = list(s = seu), selectGenesMethod = NULL,
                                 customGenelist = shared_genes,
                                 premin.spots = 1, premin.features = 1,
                                 postmin.spots = 1, postmin.features = 1)
  preobj <- AddAdjList(preobj, platform = "Visium")
  preobj <- AddParSetting(preobj, Sigma_equal = FALSE, maxIter = 30, verbose = TRUE)
  PRECAST(preobj, K = 12)
}

precast.ref <- run_precast(seu.ref)
precast.qry <- run_precast(seu.qry)

# Store PRECAST embedding in Seurat
seu.ref[["precast"]] <- CreateDimReducObject(embeddings = precast.ref@hZ$`s`, key = "PRE_", assay = DefaultAssay(seu.ref))
seu.qry[["precast"]] <- CreateDimReducObject(embeddings = precast.qry@hZ$`s`, key = "PRE_", assay = DefaultAssay(seu.qry))

# ─── Label transfer ────────────────────────────────────
seu.ref[[label_col]] <- colData(spe.ref)[[label_col]]
ref_labels <- as.character(seu.ref@meta.data[[label_col]])

anchors <- FindTransferAnchors(
  reference = seu.ref,
  query = seu.qry,
  normalization.method = "SCT",
  reduction = "precast",
  dims = 1:15
)

pred <- TransferData(
  anchorset = anchors,
  refdata = ref_labels,
  reduction = "precast",
  dims = 1:15,
  weight.reduction = seu.qry[["precast"]]
)

seu.qry <- AddMetaData(seu.qry, metadata = pred)

# ─── Push to SPE ───────────────────────────────────────
pred_lab <- seu.qry$predicted.id
pred_max <- seu.qry$prediction.score.max

colData(spe.query)[colnames(seu.qry), paste0(label_col, "_Seurat_from_", ref_id)] <- pred_lab
colData(spe.query)[colnames(seu.qry), paste0("Seurat_predScoreMax_from_", ref_id)] <- pred_max

score_cols <- grep("^prediction.score.", colnames(seu.qry@meta.data), value = TRUE)
for (sc in score_cols) {
  colData(spe)[colnames(seu.qry), paste0("Seurat_", sc, "_from_", ref_id)] <- seu.qry@meta.data[[sc]]
}

# ─── Plotting ──────────────────────────────────────────
WM_pal <- c("#D3D3D3")
green_pal <- colorRampPalette(c("#0fdb71", "#05a150"))(2)
Ce_pal <- green_pal[1]; Me_pal <- green_pal[2]
HPC_pal <- c("#FFA500")
blue_pal <- colorRampPalette(c("#03dffc", "#038cfc"))(6)
PCo_pal <- blue_pal[1]; LA_pal <- blue_pal[2]; BLVM_pal <- blue_pal[3]
BL_pal <- blue_pal[4]; BM_pal <- blue_pal[5]; BADL_pal <- blue_pal[6]
Vascular_pal <- c("#808080")
ITC_pal <- c("#ff0000ff")

pal <- c(WM_pal, Ce_pal, Me_pal, HPC_pal, PCo_pal, LA_pal, BLVM_pal, BL_pal, BM_pal, BADL_pal, Vascular_pal, ITC_pal)
names(pal) <- c("WM", "Ce", "Me", "HPC", "PCo", "LA", "BLVM", "BL", "BM", "BADL", "Vascular", "ITC")

plot_file <- here("plots", "Visium", "08_marker_genes", "label_transfer", paste0(query_id, "_label_transfer.pdf"))
dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)

pdf(plot_file, width = 10, height = 10)
p <- make_escheR(spe.query) |>
  add_fill(var = paste0(label_col, "_Seurat_from_", ref_id), point_size = 0.6) +
  scale_fill_manual(values = pal)
print(p)
dev.off()

message("Saved plot: ", plot_file)

# ─── Save ──────────────────────────────────────────────
out_file <- here("processed-data", "Visium", "08_marker_genes",
                 paste0("spe_seurat_labeltransfer_", ref_id, "_to_", query_id, ".Rdata"))
save(spe, file = out_file)

message("Saved: ", out_file)

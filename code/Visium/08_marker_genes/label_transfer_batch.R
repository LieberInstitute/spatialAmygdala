args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Please provide exactly one sample_id as argument")
query_id <- args[1]

library(Seurat)
library(SpatialExperiment)
library(Matrix)
library(scuttle)
library(here)
library(escheR)

ref_id <- "Br8325"
label_col <- "BS_manual_ITC_smoothed"

message("Running label transfer: ", ref_id, " → ", query_id)

# ─── Load data ─────────────────────────────────────────
load(here("processed-data","Visium", "06_batch_correction", "spe_harmony.Rdata"))
spe.8325 <- readRDS(here("processed-data","Visium","08_marker_genes", "spe_bs_8325_ITC.rds"))
svg_file <- read.csv(here("processed-data", "Visium", "04_feature_selection", "nnSVG_summary.csv"))
top_svg <- svg_file$gene_name[1:2000]

if (!is.null(rowData(spe)$gene_name)) rownames(spe) <- rowData(spe)$gene_name

spe.ref <- spe.8325
spe.query <- spe[, spe$sample_id == query_id]

# ─── Subset and normalize ──────────────────────────────
shared_genes <- intersect(top_svg, rownames(spe.ref))
stopifnot(length(shared_genes) >= 1000)

ref_counts <- assay(spe.ref[shared_genes, ], "counts")
qry_counts <- assay(spe.query[shared_genes, ], "counts")

ref_counts <- ref_counts[, Matrix::colSums(ref_counts) > 0]
qry_counts <- qry_counts[, Matrix::colSums(qry_counts) > 0]


options(future.globals.maxSize = 2000 * 1024^2)

seu.ref <- CreateSeuratObject(ref_counts, project = paste0("ref_", ref_id))
seu.qry <- CreateSeuratObject(qry_counts, project = paste0("qry_", query_id))
VariableFeatures(seu.ref) <- shared_genes
VariableFeatures(seu.qry) <- shared_genes

seu.ref <- SCTransform(seu.ref, verbose = FALSE, variable.features.rv.th = NULL)
seu.ref <- RunPCA(seu.ref, features = shared_genes, verbose = FALSE)

seu.qry <- SCTransform(seu.qry, verbose = FALSE, variable.features.rv.th = NULL)
seu.qry <- RunPCA(seu.qry, features = shared_genes, verbose = FALSE)

# ─── Label transfer ────────────────────────────────────
seu.ref[[label_col]] <- colData(spe.ref)[[label_col]]
ref_labels <- as.character(seu.ref@meta.data[[label_col]])

anchors <- FindTransferAnchors(
  reference = seu.ref,
  query = seu.qry,
  normalization.method = "SCT",
  dims = 1:30
)

pred <- TransferData(
  anchorset = anchors,
  refdata = ref_labels,
  dims = 1:30,
  weight.reduction = seu.qry[["pca"]]
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
pal <- c(
  WM        = "#D3D3D3",   # Light grey (WM)
  Ce        = "#1b9e77",   # Teal
  Me        = "#d95f02",   # Orange
  HPC       = "#7570b3",   # Purple
  PCo       = "#e7298a",   # Pink/magenta
  LA        = "#66a61e",   # Olive green
  BLVM      = "#e6ab02",   # Yellow-brown
  BL        = "#a6761d",   # Dark brown
  BM        = "#666666",   # Dark grey
  BADL      = "#1f78b4",   # Blue
  Vascular  = "#999999",   # Mid grey
  ITC       = "#e41a1c"    # Bright red (ITC standout)
)

plot_file <- here("plots", "Visium", "08_marker_genes", "label_transfer", paste0(query_id, "_label_transfer.pdf"))
dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)

pdf(plot_file, width = 10, height = 10)
p <- make_escheR(spe.query) |>
  add_fill(var = "BS_manual_ITC_smoothed_Seurat_from_Br8325", point_size = 1.25) +
  scale_fill_manual(values = pal)
print(p)
dev.off()

message("Saved plot: ", plot_file)

# ─── Save ──────────────────────────────────────────────
out_file <- here("processed-data", "Visium", "08_marker_genes",
                 paste0("spe_seurat_labeltransfer_", ref_id, "_to_", query_id, ".Rdata"))
save(spe, file = out_file)

message("Saved: ", out_file)

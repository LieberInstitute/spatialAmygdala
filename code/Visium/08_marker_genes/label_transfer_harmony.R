# ─── Libraries ─────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(SpatialExperiment)
  library(Matrix)
  library(scuttle)
  library(here)
  library(escheR)
})

ref_id <- "Br8325"
label_col <- "BS_manual_ITC_smoothed"

options(future.globals.maxSize = 5000 * 1024^2)

# ─── Load data ─────────────────────────────────────────
load(here("processed-data", "Visium", "06_batch_correction", "spe_harmony.Rdata"))  # loads `spe`
spe.8325 <- readRDS(here("processed-data","Visium","08_marker_genes", "spe_bs_8325_ITC.rds"))
svg_file <- read.csv(here("processed-data", "Visium", "04_feature_selection", "nnSVG_summary.csv"))
top_svg <- svg_file$gene_name[1:2000]

if (!is.null(rowData(spe)$gene_name)) rownames(spe) <- rowData(spe)$gene_name

# ─── Set up reference and query ───────────────────────
spe.ref <- spe.8325
spe.qry <- spe  # full object, includes all samples

shared_genes <- intersect(top_svg, rownames(spe))
stopifnot(length(shared_genes) >= 1000)

ref_counts <- assay(spe.ref[shared_genes, ], "counts")
qry_counts <- assay(spe.qry[shared_genes, ], "counts")
ref_counts <- ref_counts[, Matrix::colSums(ref_counts) > 0]
qry_counts <- qry_counts[, Matrix::colSums(qry_counts) > 0]

# ─── Create Seurat objects ────────────────────────────
seu.ref <- CreateSeuratObject(ref_counts, project = paste0("ref_", ref_id))
seu.qry <- CreateSeuratObject(qry_counts, project = "qry_all")
VariableFeatures(seu.ref) <- shared_genes
VariableFeatures(seu.qry) <- shared_genes

seu.ref <- SCTransform(seu.ref, ncells=15000, vars.to.regress = "nCount_RNA", verbose = FALSE, variable.features.rv.th = NULL)
seu.qry <- SCTransform(seu.qry, ncells=15000, vars.to.regress = "nCount_RNA", verbose = FALSE, variable.features.rv.th = NULL)

# ─── Assign harmonized PCA to reference ───────────────
harmony_all <- reducedDim(spe, "PCA-HARMONY_sample")
harmony_ref <- harmony_all[colnames(spe.ref), ]

seu.ref[["harmony.pca"]] <- CreateDimReducObject(
  embeddings = harmony_ref,
  assay = "SCT",
  key = "HPC_"
)
seu.ref[["pca"]] <- seu.ref[["harmony.pca"]]

# ─── Run PCA on full query (required for projection) ──
seu.qry <- RunPCA(seu.qry, features = shared_genes, verbose = FALSE)

# ─── Label transfer ───────────────────────────────────
seu.ref[[label_col]] <- colData(spe.ref)[[label_col]]
ref_labels <- as.character(seu.ref@meta.data[[label_col]])

anchors <- FindTransferAnchors(
  reference = seu.ref,
  query = seu.qry,
  normalization.method = "SCT",
  reduction = "cca",
  reference.reduction = "harmony.pca",
  dims = 1:30,
  k.anchor=50,
  k.score=50,
  max.features=500,
  mapping.score.k=TRUE
)

pred <- TransferData(
  anchorset = anchors,
  refdata = ref_labels,
  dims = 1:30,
  weight.reduction = seu.qry[["pca"]]
)
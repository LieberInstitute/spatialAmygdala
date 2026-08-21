## 02_harmony_integration.R
## Normalize the combined ITC SCE, run PCA, correct with Harmony,
## and compute UMAP on corrected embeddings.

suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
    library(scuttle)
    library(scran)
    library(scater)
    library(harmony)
    library(ggplot2)
    library(patchwork)
})

# ---- Load combined object ----
sce <- readRDS(here("processed-data", "snRNAseq", "sce_ITC_combined.rds"))

sce$ITC_type <- sce$celltype_fine

# ---- 1. Normalize ----
sce <- logNormCounts(sce)

# ---- 2. Feature selection ----
dec <- modelGeneVar(sce, block = sce$dataset)
hvgs <- getTopHVGs(dec, n = 2000)

# ---- 3. PCA ----
sce <- runPCA(sce, subset_row = hvgs, ncomponents = 30)

# ---- 4. Harmony integration ----
# Correct for both sample_id (library) and dataset
sce <- RunHarmony(
    object = sce,
    group.by.vars = c("dataset", "sample_id")
)

# ---- 5. UMAP on corrected embeddings ----
set.seed(49287)
sce <- runUMAP(sce, dimred = "HARMONY", name = "UMAP_harmony")

# Also compute uncorrected UMAP for comparison
sce <- runUMAP(sce, dimred = "PCA", name = "UMAP_uncorrected")

# ---- 6. Diagnostic plots ----
message("Generating diagnostic plots...")

# Uncorrected UMAP
p1 <- plotReducedDim(sce, "UMAP_uncorrected", colour_by = "dataset") +
    ggtitle("Uncorrected - Dataset")
p2 <- plotReducedDim(sce, "UMAP_uncorrected", colour_by = "ITC_type") +
    ggtitle("Uncorrected - Cell type")

# Corrected UMAP
p3 <- plotReducedDim(sce, "UMAP_harmony", colour_by = "dataset") +
    ggtitle("Harmony - Dataset")
p4 <- plotReducedDim(sce, "UMAP_harmony", colour_by = "ITC_type") +
    ggtitle("Harmony - Cell type")

p_combined <- (p1 | p2) / (p3 | p4)

pdf(here("plots", "snRNAseq", "ITCs", "harmony_integration_diagnostics.pdf"),
    width = 12, height = 10)
print(p_combined)

# Per-dataset facet
p5 <- plotReducedDim(sce, "UMAP_harmony", colour_by = "ITC_type") +
    facet_wrap(~sce$dataset) +
    ggtitle("Harmony UMAP by dataset")
print(p5)

# Sample-level mixing
p6 <- plotReducedDim(sce, "UMAP_harmony", colour_by = "sample_id") +
    ggtitle("Harmony - Sample") +
    theme(legend.position = "none")
print(p6)
dev.off()

message("Diagnostic plots saved.")

# ---- 7. Save ----
out_path <- here("processed-data", "snRNAseq","BICCN", "sce_ITC_harmony.rds")
saveRDS(sce, out_path)
message("Saved to: ", out_path)

sessioninfo::session_info()
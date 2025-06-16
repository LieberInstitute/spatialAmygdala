# 03_nnSVG_per-sample.R
library("SpatialExperiment")
library("here")
library("scater")
library("scran")
library("spatialLIBD")
library("nnSVG")

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
message("Running nnSVG for sample: ", sample_id)

# Save directories
plot_dir <- here("plots", "04_feature_selection")
processed_dir <- here("processed-data","Visium", "04_feature_selection")

# Load stitched object
load(here("processed-data", "Visium", "03_qc_metrics", "spe_stitched_local_outliers.Rdata"))
rownames(spe) <- rowData(spe)$gene_name

spe <- computeLibraryFactors(spe)
spe <- spe[, sizeFactors(spe) > 0]
spe <- logNormCounts(spe)

# Subset to current sample
spe_sub <- spe[, colData(spe)$sample_id == sample_id]

# Filter genes
spe_sub <- filter_genes(spe_sub, filter_genes_pcspots = 0.1)
ix_zeros <- colSums(counts(spe_sub)) == 0
if (sum(ix_zeros) > 0) {
    spe_sub <- spe_sub[, !ix_zeros]
}


set.seed(123)
spe_sub <- nnSVG(spe_sub, n_threads=10)

# Save rowData
res <- rowData(spe_sub)
saveRDS(res, file = file.path(processed_dir, paste0("nnSVG_rowdata_", sample_id, ".rds")))

# Write CSV with summary output
summary_df <- data.frame(
    gene_id = rownames(res),
    gene_name = res$gene_name,
    rank = res$rank,
    morans_I = res$morans_I,
    q_value = res$qval,
    row.names = NULL
)

csv_path <- file.path(processed_dir, paste0("nnSVG_summary_", sample_id, ".csv"))
write.csv(summary_df, csv_path, row.names = FALSE)
message("Saved summary to: ", csv_path)

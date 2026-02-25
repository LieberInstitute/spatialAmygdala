suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("RColorBrewer")
    library("ggplot2")
    library("patchwork")
    library("DeconvoBuddies")
})


output_dir <- here("plots", "Xenium", "05_spatial_clustering")

# load
spe <- readRDS(file.path("processed-data", "Xenium", "04_clustering", "Banksy", "Banksy_integrated_res2.0_collapsed_v6.rds"))
spe

# number of top genes to save
top_n <- 10

markers <- scran::findMarkers(
x         = spe,
groups    = spe$Banksy_res2.0_collapsed_v6,
test.type = "t",
pval.type = "any",    # faster, only need ranks
full.stats = FALSE,
sorted    = TRUE,
direction = "up",
assay.type = "nucleus_normcounts"
)

# get matrix of top genes per cluster, cluster by genes

top_genes <- lapply(names(markers), function(cl) {
  res <- markers[[cl]]
  # Top genes are already sorted by significance
  head(rownames(res), top_n)
})
names(top_genes) <- names(markers)

# Get unique gene set across all clusters
all_top_genes <- unique(unlist(top_genes))

# Build a matrix: clusters (rows) x genes (cols)
# Using the summary.logFC or Top rank — let's use Top rank for a clean heatmap
# Lower Top = more specific marker

marker_matrix <- matrix(NA,
  nrow = length(markers),
  ncol = length(all_top_genes),
  dimnames = list(names(markers), all_top_genes)
)

for (cl in names(markers)) {
  res <- markers[[cl]]
  shared <- intersect(all_top_genes, rownames(res))
  # Use summary.logFC (mean logFC across pairwise comparisons)
  marker_matrix[cl, shared] <- res[shared, "summary.logFC"]
}

library(pheatmap)
pdf(file.path(output_dir, "Marker_Genes_Heatmap_Banksy_res2.0_collapsed_v6.pdf"), width = 10, height = 8)
p <- pheatmap(
  marker_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "column",        # scale per gene to see relative enrichment
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  fontsize_row = 8,
  fontsize_col = 7,
  main = "Top marker genes per cluster (summary logFC)"
)
p
dev.off()


fn <- file.path(out_dir, paste0(cl, "__top", top_n, "_genes.csv"))
write.csv(wide, fn, row.names = FALSE)

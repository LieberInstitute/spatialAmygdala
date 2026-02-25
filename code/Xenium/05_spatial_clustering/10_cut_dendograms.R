#!/usr/bin/env Rscript
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(ggspavis)
  library(scCustomize)
  library(patchwork)
  library(ggplot2)
})

# ======== Parameters ========
lambda <- 0.8
k_geom <- 50
resolutions <- c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)

input_dir <- here("processed-data", "Xenium", "04_clustering", "Banksy")
coord_reset_path <- here("processed-data", "Xenium", "04_dim_reduction", "sce_combined_harmonized_singlecell.rds")
output_dir <- here("plots", "Xenium", "05_spatial_clustering")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load coordinates from original unmodified SPE
spe_orig <- readRDS(coord_reset_path)
coords_orig <- spatialCoords(spe_orig)

file_rds <- file.path(input_dir, paste0("Banksy_integrated_lambda_", lambda, "_res2_gex.rds"))

# load
spe <- readRDS(file_rds)
spe

# res 2.0
spe$custom_cluster <- factor(
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(21, 39, 4, 7), "Cluster 1",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(19,22,5,16), "Cluster 2",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(44,41,37,10,14), "Cluster 3",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(47), "Cluster 4",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(18,42,17,29,15,35), "Cluster 5",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(13,32,8,23,33,2,20,43), "Cluster 6",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(6, 1, 30), "Cluster 7",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(9, 40, 11, 12, 31), "Cluster 8",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(25), "Cluster 9",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(45, 46, 3, 36, 34, 26, 27), "Cluster 10",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(24, 28, 38), "Cluster 11",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(52, 48, 49), "Cluster 12",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(51, 54), "Cluster 13",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(53, 50, 55), "Cluster 14",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(56, 58, 57, 60), "Cluster 15",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(64, 65), "Cluster 16",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(59, 62, 61, 63), "Cluster 17",
  "Other"))))))))))))))))))

# split cluster 8
spe$custom_cluster <- factor(
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(21, 39, 4, 7), "Cluster 1",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(19,22,5,16), "Cluster 2",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(44,41,37,10,14), "Cluster 3",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(47), "Cluster 4",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(18,42,17,29,15,35), "Cluster 5",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(13,32,8,23,33,2,20,43), "Cluster 6",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(6, 1, 30), "Cluster 7",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(9), "Cluster 8",
  #ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(40), "Cluster 9",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(11), "Cluster 10",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(12, 40), "Cluster 11",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(31), "Cluster 12",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(25), "Cluster 13",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(45, 46, 3, 36, 34, 26, 27), "Cluster 14",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(24, 28, 38), "Cluster 15",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(52, 48, 49), "Cluster 16",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(51, 54), "Cluster 17",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(53, 50, 55), "Cluster 18",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(56, 58, 57, 60), "Cluster 19",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(64, 65), "Cluster 20",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2%in% c(59, 62, 61, 63), "Cluster 21",
  "Other")))))))))))))))))))))

# uncollapse Cluster 2 in res 2.0
spe$custom_cluster <- factor(
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(21, 39, 4, 7), "Cluster 1",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(19), "Cluster 2",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(22), "Cluster 3",
  #ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(5), "Cluster 4",
  #ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(16), "Cluster 5",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(44,41,37,10,14), "Cluster 6",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(47), "Cluster 7",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(18,42,17,29,15,35), "Cluster 8",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(32,8,23,33), "Cluster 10",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(2,20,43), "Cluster 11",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(6, 1, 30, 5, 16, 13), "Cluster 12",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(9), "Cluster 13",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(11), "Cluster 14",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(12, 40), "Cluster 15",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(31), "Cluster 16",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(25), "Cluster 17",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(45, 46, 3, 36, 34, 26, 27), "Cluster 18",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(24, 28, 38), "Cluster 19",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(52, 48, 49), "Cluster 20",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(51, 54), "Cluster 21",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(53, 50, 55), "Cluster 22",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(56, 58, 57, 60), "Cluster 23",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(64, 65), "Cluster 24",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(59, 62, 61, 63), "Cluster 25",
  "Other")))))))))))))))))))))))


# uncollapse Cluster 10 in res 2.0
spe$custom_cluster <- factor(
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(21, 39, 4, 7), "Cluster 1",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(19), "Cluster 2",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(22), "Cluster 3",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(44,41,37,10,14), "Cluster 4",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(47), "Cluster 5",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(18,42,17,29,15,35), "Cluster 6",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(32) , "Cluster 7",
  #ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(8), "Cluster 8",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(23), "Cluster 9",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(33), "Cluster 10",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(2,20,43,8), "Cluster 11",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(6, 1, 30, 5, 16, 13), "Cluster 12",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(9), "Cluster 13",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(11), "Cluster 14",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(12, 40), "Cluster 15",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(31), "Cluster 16",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(25), "Cluster 17",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(45, 46, 3, 36, 34, 26, 27), "Cluster 18",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(24, 28, 38), "Cluster 19",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(52, 48, 49), "Cluster 20",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(51, 54), "Cluster 21",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(53, 50, 55), "Cluster 22",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(56, 58, 57, 60), "Cluster 23",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(64, 65), "Cluster 24",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(59, 62, 61, 63), "Cluster 25",
  "Other")))))))))))))))))))))))))


  # uncollapse Cluster 18 in res 2.0
  spe$custom_cluster <- factor(
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(21, 39, 4, 7), "Cluster 1",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(19), "Cluster 2",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(22), "Cluster 3",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(44,41,37,10,14), "Cluster 4",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(47), "Cluster 5",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(18,42,17,29,15,35), "Cluster 6",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(32) , "Cluster 7",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(23), "Cluster 8",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(33), "Cluster 9",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(2,20,43,8), "Cluster 10",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(6, 1, 30, 5, 16, 13), "Cluster 11",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(9), "Cluster 12",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(11), "Cluster 13",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(12, 40), "Cluster 14",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(31), "Cluster 15",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(25), "Cluster 16",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(45), "Cluster 17",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(46), "Cluster 18",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(3), "Cluster 19",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(36), "Cluster 20",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(34), "Cluster 21",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(26), "Cluster 22",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(27), "Cluster 23",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(24, 28, 38), "Cluster 24",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(52, 48, 49), "Cluster 25",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(51, 54), "Cluster 26", 
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(53, 50, 55), "Cluster 27",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(56, 58, 57, 60), "Cluster 28",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(64, 65), "Cluster 29",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(59, 62, 61, 63), "Cluster 30",
  "Other")))))))))))))))))))))))))))))))

# collapsing similar clusters after new dendogram
  spe$custom_cluster <- factor(
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(21, 39, 4, 7, 26, 19, 22, 44,41,37,10,14), "Cluster 1",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(47), "Cluster 2",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(18,42,17,29,15,35), "Cluster 3",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(32) , "Cluster 4",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(23), "Cluster 5",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(33), "Cluster 6",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(2,20,43,8), "Cluster 7",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(6, 1, 30, 5, 16, 13), "Cluster 8",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(9), "Cluster 9",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(11), "Cluster 10",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(12, 40), "Cluster 11",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(31), "Cluster 12",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(25), "Cluster 13",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(45), "Cluster 14",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(46), "Cluster 15",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(3), "Cluster 16",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(36), "Cluster 17",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(34, 27), "Cluster 18",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(24, 28, 38), "Cluster 19",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(52, 48, 49), "Cluster 20",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(51, 54), "Cluster 21", 
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(53, 50, 55), "Cluster 22",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(56, 58, 57, 60), "Cluster 23",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(64, 65), "Cluster 24",
  ifelse(spe$clust_HARMONY_M0_lam0.8_k50_res2 %in% c(59, 62, 61, 63), "Cluster 25",
  "Other"))))))))))))))))))))))))))


num_groups <- length(unique(spe$custom_cluster))
colors <- scCustomize::scCustomize_Palette(
    num_groups = num_groups,
    ggplot_default_colors = FALSE,
    color_seed = 123
  )


# remove x_centroid and y_centroid from colData
if ("x_centroid" %in% colnames(colData(spe))) {
  colData(spe)$x_centroid <- NULL
}
if ("y_centroid" %in% colnames(colData(spe))) {
  colData(spe)$y_centroid <- NULL
}

# Plot with ggspavis
pdf_path <- file.path(output_dir, paste0("Banksy_integrated_lambda_res2.0_collapsed_v6.pdf"))
pdf(pdf_path, width = 15, height = 15)
ggspavis::plotCoords(spe,
    annotate = 'custom_cluster',
    in_tissue = NULL,
    sample_id = "brnum"
  ) +
    scale_color_manual(values = colors) +
    ggtitle(paste0("res=2.0 collapsed")) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    guides(color = guide_legend(nrow = 2, 
                    byrow = TRUE, 
                    override.aes = list(size=5)
                    )
            )
dev.off()




# ======= Dendrograms of new clusters =======

library(dreamlet)

# drop polygons colData if exists
if ("polygons" %in% colnames(colData(spe))) {
  colData(spe)$polygons <- NULL
}

# aggregate by cluster + sample
pb <- aggregateToPseudoBulk(
  spe,
  assay = "counts",
  cluster_id = "custom_cluster",
  sample_id = "sample_id",
  verbose = FALSE
)

# build dendrogram
hcl <- buildClusterTreeFromPB(pb, method = "ward.D")

# save plot
pdf_path <- file.path(output_dir, paste0("Dendrogram_res2.0_collapsed_v5.pdf"))
pdf(pdf_path, width = 8, height = 4)
plot(hcl, hang = -1, main = paste0("res = 2.0 collapsed"))
dev.off()



# ===== Violins and bar plots of library size and n_cells per cluster =====

pdf_path <- file.path(output_dir, paste0("Library_size_per_cluster_res2.0_collapsed_v5.pdf"))
pdf(pdf_path, width = 8, height = 4)
scater::plotColData(spe, x="custom_cluster", y="total", color_by="custom_cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Library size per cluster") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(trans = "log10")+
  scale_x_discrete(limits = hcl$labels[hcl$order])
dev.off()



# save custom clusters as Banksy_res2.0_collapsed_v6 in colData
colData(spe)$Banksy_res2.0_collapsed_v6 <- spe$custom_cluster
spe$custom_cluster <- NULL

# save spe
saveRDS(spe, file.path("processed-data", "Xenium", "04_clustering", "Banksy", "Banksy_integrated_res2.0_collapsed_v6.rds"))
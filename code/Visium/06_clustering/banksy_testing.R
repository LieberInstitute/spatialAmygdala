library("here")
library("SpatialExperiment")
library("SummarizedExperiment")
library("SingleCellExperiment")
library("spatialLIBD")
library("ggplot2")
library("patchwork")
library("Banksy")
library("Seurat")
library("scater")
library("scran")
library("harmony")
library("escheR")


# Save directories
plot_dir = here("plots", "08_clustering", "banksy")
processed_dir = here("processed-data","08_clustering")

load(here("processed-data","Visium", "06_dim_reduction", "spe_stitched_pca.Rdata"), verbose = TRUE)
spe
# class: SpatialExperiment 
# dim: 25801 308044 
# metadata(0):
# assays(2): counts logcounts
# rownames(25801): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(6): source type ... gene_name Symbol.uniq
# colnames(308044): AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 ...
#   TTGTTTCCATACAACT-1 TTGTTTGTGTAAATTC-1
# colData names(36): sample_id in_tissue ... local_outliers sizeFactor
# reducedDimNames(4): 10x_pca 10x_tsne 10x_umap PCA
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

# subset to fourth donor (Br9280)
spe.subset <- spe[, spe$brnum == "Br8325"]

colnames(colData(spe.subset))
#unique(spe$exclude_overlapping)
#[1] FALSE  TRUE    NA

# eclude overlapping and NAs
spe.subset <- spe.subset[, !spe.subset$exclude_overlapping & !is.na(spe.subset$exclude_overlapping)]

# renormalization
spe.subset <- logNormCounts(spe.subset)

# get HVGs
dec <- modelGeneVar(spe.subset)
chosen <- getTopHVGs(dec, n=4000)

# subset to hvgs
spe.subset <- spe.subset[chosen, ]


# == Banksy ==sq
lambda <- 0.4
k_geom <- 18
npcs <- 20
aname <- "logcounts"

spe.subset <- Banksy::computeBanksy(spe.subset, assay_name = aname, k_geom = k_geom)

set.seed(1000)
spe.subset <- Banksy::runBanksyPCA(spe.subset, lambda = lambda, npcs = npcs, group = "capture_area")

# # drop PCA, rename PCA_M0_lam0.4 to PCA
# reducedDim(spe.subset, "PCA") <- NULL
# reducedDim(spe.subset, "HARMONY") <- NULL
# reducedDim(spe.subset, "PCA") <- reducedDim(spe.subset, "PCA_M0_lam0.4")
# reducedDim(spe.subset, "PCA_M0_lam0.4") <- NULL

# spe.subset <- RunHarmony(spe.subset, "capture_area")

# set.seed(1000)
# spe.subset$clust_HARMONY_k50_res0.8 <- NULL
# spe.subset <- Banksy::clusterBanksy(spe.subset, lambda = lambda, npcs = npcs, resolution = 0.4, dimred = "HARMONY")



# drop PCA, rename PCA_M0_lam0.4 to PCA

set.seed(1000)
#spe.subset$clust_HARMONY_k50_res0.8 <- NULL
spe.subset <- Banksy::clusterBanksy(spe.subset, lambda = lambda, npcs = npcs, resolution = 0.4, dimred = "PCA_M0_lam0.4")

# drop duplicate columns
colData(spe.subset) <- colData(spe.subset)[, !duplicated(colnames(colData(spe.subset)))]


pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(spe.subset$clust_PCA_M0_lam0.4_k50_res0.4)))
png(file = here("plots", "Visium", "08_clustering", "banksy", "Br8325_Banksy_lambda_0.4_susbset_multisample_k18_res0.4.png"), width=10, height=10, units="in", res=300)
make_escheR(spe.subset) |>
    add_fill(var="clust_PCA_M0_lam0.4_k50_res0.4") +
    scale_fill_manual(values=pal) 
dev.off()




saveRDS(spe.subset, here("processed-data", "Visium","08_clustering", "Banksy", "Br8325_Banksy_lambda_0.8_stitched.rds"))














# # ===== Create sample list of SPEs =====
# sample_names <- unique(spe$sample_id)
# spe_list <- lapply(sample_names, function(x) spe[, spe$sample_id == x])


# # ===== Preprocessing =====
# #' Normalize data
# spe_list <- lapply(spe_list, function(x) {
#     logNormCounts(x)
# })

# #' Compute HVGs
# var <- lapply(spe_list, function(x) {
#      modelGeneVar(x)
# })

# # get top HVGs
# hvg <- lapply(var, function(x) {
#     getTopHVGs(x, n=500)
# })
# hvg <- Reduce(union, hvg)


# # Subset each SingleCellExperiment to HVGs
# spe_list_hvg <- lapply(seq_along(spe_list), function(i) {
#     spe <- spe_list[[i]]
#     spe[hvg, ]
# })

# # ===== Run Banksy =====
# compute_agf <- TRUE
# k_geom <- 6 # first-order neighbors
# spe_list <- lapply(spe_list_hvg, computeBanksy, assay_name = "logcounts", 
#                    compute_agf = compute_agf, k_geom = k_geom)

# # merge samples to perform joint dimensional reduction and clustering
# spe_joint <- do.call(cbind, spe_list)
# rm(spe_list)
# invisible(gc())

# # lamba function determines neighborhood weighting. groups will be z-scaled seperately. 
# lambda <- 0.4
# use_agf <- TRUE
# spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda, group = "brnum", seed = 1000)

# set.seed(1000)
# harmony_embedding <- RunHarmony(
#     data_mat = reducedDim(spe_joint, "PCA_M1_lam0.4"),
#     meta_data = colData(spe_joint),
#     vars_use = c("sample_id", "brnum"),
#     do_pca = FALSE,
#     max_iter = 20,
#     verbose = FALSE
# )
# reducedDim(spe_joint, "Harmony_BANKSY") <- harmony_embedding

# # UMAP on banksy embedding
# spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
# spe_joint <- runBanksyUMAP(spe_joint, dimred = "Harmony_BANKSY")

# reducedDims(spe_joint)
# # 10x_pca 10x_tsne 10x_umap PCA PCA_M0_lam0.2 UMAP_M0_lam0.2

# png(here(plot_dir, "banksy_joint_umap.png"), width = 20, height = 10, units = "in", res = 300)
# p1 <- plotReducedDim(spe_joint, dimred="UMAP_M1_lam0.4", colour_by="brnum")
# p2 <- plotReducedDim(spe_joint, dimred="UMAP_Harmony_BANKSY", colour_by="brnum")
# p1+p2
# dev.off()

# png(here(plot_dir, "banksy_joint_umap_broad_markers.png"), width = 20, height = 10, units = "in", res = 300)
# p1 <- plotReducedDim(spe_joint, dimred="UMAP_Harmony_BANKSY", colour_by="SNAP25")
# p2 <- plotReducedDim(spe_joint, dimred="UMAP_Harmony_BANKSY", colour_by="SLC17A7")
# p3 <- plotReducedDim(spe_joint, dimred="UMAP_Harmony_BANKSY", colour_by="GAD1")

# p4 <- plotReducedDim(spe_joint, dimred="UMAP_Harmony_BANKSY", colour_by="sum_umi")
# p5 <- plotReducedDim(spe_joint, dimred="UMAP_Harmony_BANKSY", colour_by="sum_gene")
# p6 <- plotReducedDim(spe_joint, dimred="UMAP_Harmony_BANKSY", colour_by="expr_chrM_ratio")
# (p1+p2+p3)/ (p4+p5+p6) 
# dev.off()

# # obtain cluster labels across all samples
# res <- 0.7
# spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)
# cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
# #spe_joint <- connectClusters(spe_joint)

# # split sample into their own SPE again
# spe_list <- lapply(sample_names, function(x) spe_joint[, spe_joint$sample_id == x])
# rm(spe_joint)
# invisible(gc())

# # Optionally, smooth the cluster labels 
# spe_list <- lapply(spe_list, smoothLabels, cluster_names = cnm, k = 6L, verbose = FALSE)
# names(spe_list) <- paste0("sample_", sample_names)

# spe_joint <- do.call(cbind, spe_list)

# # drop dup columns
# colData(spe_joint) <- colData(spe_joint)[, !duplicated(colnames(colData(spe_joint)))]


# # make distinct color palette for each cluster using Rcolorbewer
# #subset to just just one brnum 8325
# spe.8325 <- spe_joint[, spe_joint$brnum == "Br8325"]

# cols <- Polychrome::palette36.colors(length(unique(spe.8325$clust_M1_lam0.4_k50_res0.7)))
# names(cols) <- 1:length(cols)
# cols

# p_list <- vis_grid_clus(
#     spe.8325,
#     clustervar= "clust_M1_lam0.4_k50_res0.7",
#     spatial = FALSE,
#     auto_crop = TRUE,
#     return_plots = TRUE,
#     pdf_file = NULL,
#     colors = cols
# )
# plot_list_reordered <- c(p_list[4], p_list[6], p_list[5],
#                          p_list[3], p_list[8], p_list[7],
#                          p_list[2], p_list[1])
# cowplot::plot_grid(plotlist = plot_list_reordered, ncol = 3)
# ggsave(here(plot_dir, paste0("Br8325_banksy_grid_spotplot.png")), width = 20, height = 20)

# # save spe
# save(spe_joint, file = here(processed_dir, "spe_banksy.Rdata"))
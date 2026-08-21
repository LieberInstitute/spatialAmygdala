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

out_dir <- here("processed-data", "Visium", "08_marker_genes")

# load
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
spe 
# class: SpatialExperiment 
# dim: 36601 226181 
# metadata(0):
# assays(2): counts logcounts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type gene_search
# colnames(226181): AAACAAGTATCTCCCA-1_V13Y24-346_A1
#   AAACAATCTACTAGCA-1_V13Y24-346_A1 ... TTGTTTCATTAGTCTA-1_V13B23-407_C1
#   TTGTTTCCATACAACT-1_V13B23-407_C1
# colData names(59): sample_id in_tissue ... ITC_gmm2_smooth
#   BS_k16_relabel_ITC_smooth
# reducedDimNames(2): PCA PCA-HARMONY_sample
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


# ===== Pseudobulk marker genes using spatialLIBD ======


# renormalize
markers <- scran::findMarkers(
  spe,
  groups=spe$BS_k16_Semisupervised_wAI,
  test.type = c("wilcox"),
  pval.type = c("all"),
  full.stats = TRUE,
  sorted = TRUE,
  direction="up"
)

# get top 10 gene names for each cluster
top_markers <- lapply(markers, function(x) {
  top10 <- x[1:10,]
  top10_genes <- rownames(top10)
  top10_genes
})
# $`1`
#  [1] "PENK"    "SYNPR"   "TMEM272" "GAD2"    "SLC32A1" "SLC35F1" "GPR88"  
#  [8] "NLRP1"   "TSHZ1"   "GAD1"   

# save markers as csv
saveRDS(markers, file=here("processed-data","Visium","08_marker_genes", "markers_bs_final_ITC_smoothed_wBLVM.rds"))


# csv of top 100 markers per cluster

top100 <- lapply(markers, function(x) {
  top100 <- x[1:100,]
  top100_genes <- rownames(top100)
  top100_genes
})

# convert to data frame
top100_df <- do.call(rbind, lapply(names(top100), function(cluster) {
  data.frame(
    cluster = cluster,
    gene = top100[[cluster]],
    stringsAsFactors = FALSE
  )
}))

# write to csv  
write.csv(top100_df, file=here("processed-data","Visium","08_marker_genes", "top100_markers_bs_final_ITC_smoothed_wBLVM.csv"), row.names = FALSE)



# ======== Heatmap plotting ========
# blue to red coloscale
# pdf(here("plots", "Visium", "08_marker_genes","final_clusters", "amygdala_heatmap_known_markers.pdf"), width = 4, height = 3)
# p <- scater::plotGroupedHeatmap(spe, known, group="BS_k16_relabel_ITC_smooth", center=TRUE, scale=TRUE, treeheight_row = 0, treeheight_col = 0)
# print(p)
# dev.off()


# plot grouped heatmap of the top 5 markers per cluster
top5_genes <- unique(top100_df$gene[ave(as.numeric(top100_df$cluster), top100_df$cluster, FUN = function(x) seq_along(x) <=10) == 1])
pdf(here("plots", "Visium", "08_marker_genes","final_clusters" ,"grouped_heatmap_top5markers_wBLVM.pdf"), width = 6, height = 15)
p <- scater::plotGroupedHeatmap(spe, top5_genes, group="BS_k16_Semisupervised_wAI", center=TRUE, scale=TRUE, treeheight_row = 0, treeheight_col = 0)
print(p)
dev.off()
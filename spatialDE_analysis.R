setwd('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("BayesSpace")
    library("RColorBrewer")
    library("ggplot2")
    library("gridExtra")
    library("patchwork")
})

# Save directories
plot_dir = here("plots", "09_xenium_panel")
processed_dir = here("processed-data","09_xenium_panel")

load(here("processed-data","08_clustering", "BayesSpace", "spe_clusters_k10.Rdata"), verbose = TRUE)
spe
# class: SpatialExperiment 
# dim: 28412 29885 
# metadata(2): BayesSpace.data chain.h5
# assays(2): counts logcounts
# rownames(28412): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(7): source type ... Symbol.uniq is.HVG
# colnames(29885): AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 ...
# TTGTTTCATTAGTCTA-1 TTGTTTCCATACAACT-1
# colData names(47): sample_id in_tissue ... cluster.init spatial.cluster
# reducedDimNames(4): 10x_pca 10x_tsne 10x_umap PCA
# mainExpName: NULL
# altExpNames(0):
#     spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


# Drop any duplicate coolData
colData(spe) <- colData(spe)[ , !duplicated(colnames(colData(spe)))]
colnames(colData(spe))



# ======= Visualize spatial domains =======

p_list <- vis_grid_clus(
    spe,
    clustervar= "spatial.cluster",
    spatial = FALSE,
    auto_crop = TRUE,
    return_plots = TRUE,
    pdf_file = NULL,
)
plot_list_reordered <- c(p_list[5], p_list[6], p_list[4],
                         p_list[7], p_list[8], p_list[3],
                         p_list[1], p_list[2])
cowplot::plot_grid(plotlist = plot_list_reordered, ncol = 3)
ggsave(here(plot_dir, "BayesSpace_k10_domains.pdf"), width = 20, height = 20)


# ======= Annotate spatial domains =======

# relabel spe$spatial.cluster
spe$spatial.cluster <- factor(spe$spatial.cluster)
levels(spe$spatial.cluster) <- c('WM.1', 'aBA', 'LA.1', 'BA', 'vmBA', 'GABA.n', 'LA.2', 'EC', 'WM.2', 'Endo')
unique(spe$spatial.cluster)

p_list <- vis_grid_clus(
    spe,
    clustervar= "spatial.cluster",
    spatial = FALSE,
    auto_crop = TRUE,
    return_plots = TRUE,
    pdf_file = NULL,
)
plot_list_reordered <- c(p_list[5], p_list[6], p_list[4],
                         p_list[7], p_list[8], p_list[3],
                         p_list[1], p_list[2])
cowplot::plot_grid(plotlist = plot_list_reordered, ncol = 3)
ggsave(here(plot_dir, "BayesSpace_k10_domains_annotated.pdf"), width = 20, height = 20)





features = c("PDYN", "CYP26B1", "RGS4", "LAMP5") #, "COL25A1", "GULP1")

p <- scater::plotExpression(spe, features, x='spatial.cluster', colour_by='spatial.cluster', show_violin=TRUE, show_median=TRUE, ncol=1)

pdf(width=7, height=7, here(plot_dir, "BLA_subset_expressionPlots_top_marker_k10.pdf"))
p
dev.off()



# ======= Marker gene analysis =======
library(scran)

# find and remove MT genes
mt <- grep("^MT-", rownames(spe))
spe <- spe[-mt,]


# find marker genes
marker.info <- scoreMarkers(spe, spe$spatial.cluster)
marker.info

# ======= Visualize marker genes =======

# Note: There are many different metrics given to us by marker.info. Some are more sensitive to the magnitude
# of expression than others. For our purposes, we likely want genes that have decent overall express. 
# let's plot a few different metrics to see which works best

cluster_id <- levels(unique(spe$spatial.cluster))

# loop through cluster ids
for (i in 1:length(cluster_id)) {
    
    # just look at top 5 genes for one spatial cluster as a sanity check
    chosen <- marker.info[[cluster_id[i]]]

    plot_dir = here("plots", "09_xenium_panel", "amygdala_subregions", cluster_id[i])
    
    #make sure folder exists, if not, make it
    if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
    }
    
    # output for keeping track with loop
    print(paste("Now printing marker genes for:", cluster_id[i]))
    
    # mean AUC
    mean_AUC <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
    head(mean_AUC[,1:4]) # showing basic stats only, for brevity.
    
    # median AUC
    median_AUC <- chosen[order(chosen$median.AUC, decreasing=TRUE),]
    head(median_AUC[,1:4]) # showing basic stats only, for brevity.
    
    # mean Cohen d
    mean_Cohen_d <- chosen[order(chosen$mean.logFC.cohen, decreasing=TRUE),]
    head(mean_Cohen_d[,1:4]) # showing basic stats only, for brevity.
    
    # median Cohen d
    median_Cohen_d <- chosen[order(chosen$median.logFC.cohen, decreasing=TRUE),]
    head(median_Cohen_d[,1:4]) # showing basic stats only, for brevity.
    
    # min Cohen d
    min_Cohen_d <- chosen[order(chosen$min.logFC.cohen, decreasing=TRUE),]
    head(min_Cohen_d[,1:4]) # showing basic stats only, for brevity.
    
    # min AUC
    min_AUC <- chosen[order(chosen$min.AUC, decreasing=TRUE),]
    head(min_AUC[,1:4]) # showing basic stats only, for brevity.
    
    # visualize mean AUC
    features <- c(rownames(mean_AUC[1:5,]))
    p <- plotExpression(spe, features=features, x='spatial.cluster', colour_by='spatial.cluster', show_violin=TRUE, show_median=TRUE, ncol=1)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_mean_AUC.pdf"))
    print(p)
    dev.off()
    
    # visualize median AUC
    features <- c(rownames(median_AUC[1:5,]))
    p <- plotExpression(spe, features=features, x='spatial.cluster', colour_by='spatial.cluster', show_violin=TRUE, show_median=TRUE, ncol=1)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_median_AUC.pdf"))
    print(p)
    dev.off()
    
    # visualize mean Cohen d
    features <- c(rownames(mean_Cohen_d[1:5,]))
    p <- plotExpression(spe, features=features, x='spatial.cluster', colour_by='spatial.cluster', show_violin=TRUE, show_median=TRUE, ncol=1)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_aBA_mean_Cohen_d.pdf"))
    print(p)
    dev.off()
    
    # visualize median Cohen d
    features <- c(rownames(median_Cohen_d[1:5,]))
    p <- plotExpression(spe, features=features, x='spatial.cluster', colour_by='spatial.cluster', show_violin=TRUE, show_median=TRUE, ncol=1)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_aBA_median_Cohen_d.pdf"))
    print(p)
    dev.off()
    
    # visualize min Cohen d
    features <- c(rownames(min_Cohen_d[1:5,]))
    p <- plotExpression(spe, features=features, x='spatial.cluster', colour_by='spatial.cluster', show_violin=TRUE, show_median=TRUE, ncol=1)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_aBA_min_Cohen_d.pdf"))
    print(p)
    dev.off()
    
    # visualize min AUC
    features <- c(rownames(min_AUC[1:5,]))
    p <- plotExpression(spe, features=features, x='spatial.cluster', colour_by='spatial.cluster', show_violin=TRUE, show_median=TRUE, ncol=1)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_aBA_min_AUC.pdf"))
    print(p)
    dev.off()
}
    

    

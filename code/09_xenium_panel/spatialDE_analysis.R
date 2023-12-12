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
plot_dir = here("plots", "09_xenium_panel", "amygdala_subregions")
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


# subset spe to just spatial clusters that contain the string LA or BL
spe.subset <- spe[,grepl("LA|BA", spe$spatial.cluster)]
unique(spe.subset$spatial.cluster)
# [1] vmBA   LA.2   LA.1   aBA    GABA.n BA    
# Levels: WM.1 aBA LA.1 BA vmBA GABA.n LA.2 EC WM.2 Endo

#reset levels of spatial clusters
spe.subset$spatial.cluster <- factor(spe.subset$spatial.cluster, levels = c("LA.1", "LA.2", "aBA", "BA", "vmBA", "GABA.n"))

# collapse LA.1 and LA.2 into one cluster by changing labels to LA
spe.subset$spatial.cluster <- plyr::mapvalues(spe.subset$spatial.cluster, from = c("LA.1", "LA.2"), to = c("LA", "LA"))

# find marker genes
marker.info <- findMarkers(spe.subset, spe.subset$spatial.cluster, test="binom", direction="up", lfc=2)
marker.info

# get cluster id for loop
cluster_id <- unique(levels(spe.subset$spatial.cluster))

# loop through cluster ids
for (i in 1:length(cluster_id)) {
    
    # just look at top 5 genes for one spatial cluster as a sanity check
    chosen <- marker.info[[i]]

    plot_dir = here("plots", "09_xenium_panel", "amygdala_subregions", cluster_id[i])
    
    #make sure folder exists, if not, make it
    if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
    }
    
    # output for keeping track with loop
    print(paste("Now printing marker genes for:", cluster_id[i]))
    
    
    # visualize mean AUC
    top10 <- c(rownames(chosen[1:10,]))
    p <- scater::plotExpression(spe, features=top10, x='spatial.cluster', colour_by='spatial.cluster', show_violin=TRUE, show_median=TRUE, ncol=2)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_Top10.pdf"))
    print(p)
    dev.off()
    
}
    

# ========= Volano plots ==========
library(EnhancedVolcano)

# find marker genes
marker.info <- findMarkers(spe.subset, spe.subset$spatial.cluster, test="binom", direction="up")
marker.info

# Loop through each unique cluster ID
for(cluster_id in unique(spe.subset$spatial.cluster)){
    
    # set plot dir
    plot_dir = here("plots", "09_xenium_panel", "amygdala_subregions", cluster_id)
    
    # Filter the all_markers data frame for the current cluster
    cluster_markers <- as.data.frame(marker.info[[cluster_id]])
    
    volcano_plot <- EnhancedVolcano(
        cluster_markers,
        lab = row.names(cluster_markers),
        x = 'summary.logFC',
        y = 'p.value',
        title = paste('Volcano plot of Cluster', cluster_id, 'markers'),
        FCcutoff = 2,
        pCutoff = 10e-50,
        drawConnectors = TRUE,
        widthConnectors = 0.75
    )
    
    
    # Optionally, save the volcano plot to a file
    ggsave(filename = here(plot_dir, paste0("volcano_plot_cluster_", cluster_id, ".png")), plot = volcano_plot)
    
}


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
#mt <- grep("^MT-", rownames(spe))
#spe <- spe[!mt,]


# subset spe to just spatial clusters that contain the string LA or BL
spe.subset <- spe[,grepl("LA|BA", spe$spatial.cluster)]
unique(spe.subset$spatial.cluster)
# [1] vmBA   LA.2   LA.1   aBA    GABA.n BA    
# Levels: WM.1 aBA LA.1 BA vmBA GABA.n LA.2 EC WM.2 Endo

# also drop the GABA.n cluster
spe.subset <- spe.subset[,spe.subset$spatial.cluster != "GABA.n"]

#reset levels of spatial clusters
spe.subset$spatial.cluster <- factor(spe.subset$spatial.cluster, levels = c("LA.1", "LA.2", "aBA", "BA", "vmBA"))

# collapse LA.1 and LA.2 into one cluster by changing labels to LA
spe.subset$spatial.cluster <- plyr::mapvalues(spe.subset$spatial.cluster, from = c("LA.1", "LA.2"), to = c("LA", "LA"))

# find marker genes
marker.info <- findMarkers(spe.subset, spe.subset$spatial.cluster, test="wilcox", direction="up", lfc=, add.summary=TRUE)
marker.info

# for each region, get only markers with p.value < .05 & self.average > 1
for (i in 1:length(marker.info)) {
    marker.info[[i]] <- marker.info[[i]][marker.info[[i]]$Top < 50 & marker.info[[i]]$self.average > 2.5,]
}

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




# ======== Getting nnSVGs to round out the final list =========

# read in current list
custom_markers <- read.csv(here(processed_dir, "Amygdala_Xenium_panel_3.0.csv"), header = TRUE, stringsAsFactors = FALSE)
head(custom_markers)
# # A tibble: 6 × 2
# Region Gene   
# <chr>  <chr>  
#     1 LA     SLC30A3
# 2 LA     CYP26B1
# 3 LA     NCALD  
# 4 LA     NEUROD6
# 5 LA     CPNE4  
# 6 BA     KLK7   

# read in nnSVGs
SVGs <- read.csv(here("processed-data","05_feature_selection", "nnSVG_summary.csv"), header = TRUE, stringsAsFactors = FALSE)
#   gene_id gene_name gene_type overall_rank average_rank n_withinTop100
# 1     MBP       MBP      gene          1.5        2.000              8
# 2     AVP       AVP      gene          1.5        2.000              1
# 3    NGFR      NGFR      gene          3.0        5.000              1
# 4    ENC1      ENC1      gene          4.0        5.750              8
# 5  SNAP25    SNAP25      gene          5.0        7.625              8

# read in base panel 
base_panel <- read.csv(here(processed_dir, "Xenium_hBrain_v1_metadata.csv"), header = TRUE, stringsAsFactors = FALSE)
head(base_panel)
# Genes      Ensembl_ID Num_Probesets Codewords Annotation
# 1    ABCC9 ENSG00000069431             8         1       VLMC
# 2 ADAMTS12 ENSG00000151388             8         1       VLMC
# 3 ADAMTS16 ENSG00000145536             8         1      L4 IT
# 4  ADAMTS3 ENSG00000156140             8         1      L6 IT
# 5   ADRA1A ENSG00000120907             8         1       Sncg
# 6   ADRA1B ENSG00000170214             8         1       Sncg

# === Getting new genes ===

# get top 50 SVGs that are not in custom_markers or base_panel
new_genes <- SVGs[!SVGs$gene_name %in% custom_markers$Gene & !SVGs$gene_name %in% base_panel$Genes,]
top_genes <- new_genes[1:100,]
head(top_genes)


# remove any custom_markers that are in base_panel
custom_markers <- custom_markers[!custom_markers$Gene %in% base_panel$Genes,]


# add gene Ensemble_ID column (rowData(spe)$gene_id) for all custom_markers
custom_markers$Ensembl_ID <- rowData(spe)$gene_id[match(custom_markers$Gene, rowData(spe)$gene_name)]

# change Region column name to Annotation
colnames(custom_markers)[1] <- "Annotation"
custom_markers

# save as final xenium probe lsit
write.csv(custom_markers, here(processed_dir, "Amygdala_Xenium_panel_4.0_final.csv"), row.names = FALSE)

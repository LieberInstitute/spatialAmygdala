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
plot_dir = here("plots", "09_xenium_panel", "snRNA-seq")
processed_dir = here("processed-data","09_xenium_panel")

load(here("processed-data","snRNA-seq", "sce_annotated.rda"), verbose = TRUE)
sce
# class: SingleCellExperiment 
# dim: 36601 19682 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(36601): MALAT1 ERBB4 ... AC136616.2 AC141272.1
# rowData names(8): source type ... Symbol.uniq binomial_deviance
# colnames(19682): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
# colData names(24): Sample Barcode ... celltype annotation
# reducedDimNames(4): PCA TSNE UMAP HARMONY
# mainExpName: NULL
# altExpNames(0):


# Drop any duplicate coolData
colData(sce) <- colData(sce)[ , !duplicated(colnames(colData(sce)))]
colnames(colData(sce))
# [1] "Sample"                "Barcode"               "sum"                  
# [4] "detected"              "subsets_Mito_sum"      "subsets_Mito_detected"
# [7] "subsets_Mito_percent"  "total"                 "high_mito"            
# [10] "low_lib"               "low_genes"             "discard_auto"         
# [13] "doubletScore"          "sizeFactor"            "k_25_label"           
# [16] "k_50_label"            "celltype"              "annotation"    


unique(sce$annotation)
# [1] "Inh_LAMP5_2" "Exc_09"      "Inh_3"       "Exc_03"      "Inh_LAMP5_1"
# [6] "Inh_ITC_1"   "Exc_01"      "Inh_SST"     "Exc_05"      "Inh_PVALB"  
# [11] "OPC"         "Exc_04"      "Oligo"       "Micro"       "Inh_2"      
# [16] "Astro_5"     "Exc_12"      "Inh_5"       "Inh_VIP_1"   "Inh_4"      
# [21] "Astro_3"     "Inh_VIP_2"   "Inh_1"       "Endo"        "Inh_ITC_2"  
# [26] "Exc_08"      "Exc_07"      "Exc_11"      "Astro_1"     "Exc_02"     
# [31] "Exc_06"      "Exc_10"      "Astro_2"     "Astro_4"  

unique(sce$celltype)
# [1] "Inhib"        "Excit"        "Non-neuronal"



# ======= Marker gene analysis =======
library(scran)

# find and remove MT genes
mt <- grep("^MT-", rownames(sce))
sce <- sce[-mt,]


# subset sce to just spatial clusters that contain the string LA or BL
sce.subset <- sce[,grepl("Excit", sce$celltype)]
unique(sce.subset$celltype)
# [1] "Excit"

unique(sce.subset$annotation)
# [1] "Exc_09" "Exc_03" "Exc_01" "Exc_05" "Exc_04" "Exc_12" "Exc_08" "Exc_07"
# [9] "Exc_11" "Exc_02" "Exc_06" "Exc_10"

# find marker genes
marker.info <- findMarkers(sce.subset, sce.subset$annotation, test="binom", direction="up", lfc=2)
marker.info

# get cluster id for loop
cluster_id <- unique(levels(as.factor(sce.subset$annotation)))

# loop through cluster ids
for (i in 1:length(cluster_id)) {
    
    # just look at top 5 genes for one spatial cluster as a sanity check
    chosen <- marker.info[[i]]

    plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", cluster_id[i])
    
    #make sure folder exists, if not, make it
    if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
    }
    
    # output for keeping track with loop
    print(paste("Now printing marker genes for:", cluster_id[i]))
    
    
    # visualize mean AUC
    top10 <- c(rownames(chosen[1:10,]))
    p <- scater::plotExpression(sce.subset, features=top10, x='annotation', colour_by='annotation', show_violin=TRUE, show_median=TRUE, ncol=2)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_Top10.pdf"))
    print(p)
    dev.off()
    
}
    


# =========== Intercalated cell clusters ===========

# subset sce to just spatial clusters that contain the string ITC
sce.subset <- sce[,grepl("Inhib", sce$celltype)]

unique(sce.subset$celltype)
# [1] "Inhib"

unique(sce.subset$annotation)
# [1] "Inh_LAMP5_2" "Inh_3"       "Inh_LAMP5_1" "Inh_ITC_1"   "Inh_SST"    
# [6] "Inh_PVALB"   "Inh_2"       "Inh_5"       "Inh_VIP_1"   "Inh_4"      
# [11] "Inh_VIP_2"   "Inh_1"       "Inh_ITC_2" 

# find marker genes
marker.info <- findMarkers(sce.subset, sce.subset$annotation, test="binom", direction="up", lfc=2)
marker.info

# get cluster id for loop
cluster_id <- unique(levels(as.factor(sce.subset$annotation)))

# loop through cluster ids
for (i in 1:length(cluster_id)) {
    
    # just look at top 5 genes for one spatial cluster as a sanity check
    chosen <- marker.info[[i]]

    plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", cluster_id[i])
    
    #make sure folder exists, if not, make it
    if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
    }
    
    # output for keeping track with loop
    print(paste("Now printing marker genes for:", cluster_id[i]))
    
    
    # visualize mean AUC
    top10 <- c(rownames(chosen[1:10,]))
    p <- scater::plotExpression(sce.subset, features=top10, x='annotation', colour_by='annotation', show_violin=TRUE, show_median=TRUE, ncol=2)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_Top10.pdf"))
    print(p)
    dev.off()
    
}


# =====+ ITC 1 vs ITC 2 ======


# run find markers just comparing ITC 1 and 2
sce.subset <- sce[,grepl("ITC", sce$annotation)]

unique(sce.subset$celltype)
# [1] "Inhib"

unique(sce.subset$annotation)
# [1] "Inh_ITC_1" "Inh_ITC_2"

# find marker genes
marker.info <- marker.info <- findMarkers(sce.subset, sce.subset$annotation, test="t")
marker.info

# get cluster id for loop
cluster_id <- unique(levels(as.factor(sce.subset$annotation)))

# loop through cluster ids
for (i in 1:length(cluster_id)) {
    
    # just look at top 5 genes for one spatial cluster as a sanity check
    chosen <- marker.info[[i]]

    plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", "1v1_ITC", cluster_id[i])
    
    #make sure folder exists, if not, make it
    if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
    }
    
    # output for keeping track with loop
    print(paste("Now printing marker genes for:", cluster_id[i]))
    
    
    # visualize mean AUC
    top10 <- c(rownames(chosen[1:10,]))
    p <- scater::plotExpression(sce.subset, features=top10, x='annotation', colour_by='annotation', show_violin=TRUE, show_median=TRUE, ncol=2)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_Top10.pdf"))
    print(p)
    dev.off()
    
}




# === Volano plots ===
library(EnhancedVolcano)

marker.info <- findMarkers(sce.subset, sce.subset$annotation, test="t")

# Loop through each unique cluster ID
for(cluster_id in unique(sce.subset$annotation)){

    # set plot dir
    plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", "1v1_ITC", cluster_id)

    # Filter the all_markers data frame for the current cluster
    cluster_markers <- as.data.frame(marker.info[[cluster_id]])

    volcano_plot <- EnhancedVolcano(
        cluster_markers,
        lab = row.names(cluster_markers),
        x = 'summary.logFC',
        y = 'p.value',
        title = paste('Volcano plot of Cluster', cluster_id, 'markers'),
        FCcutoff = 2,
        pCutoff = 10e-2,
        drawConnectors = TRUE,
        widthConnectors = 0.75
    )


    # Optionally, save the volcano plot to a file
    ggsave(filename = here(plot_dir, paste0("volcano_plot_cluster_", cluster_id, ".png")), plot = volcano_plot)

}


# ========== LAMP5 1 va LAMP5 2 ===========

# run find markers just comparing LAMP5 1 and 2
sce.subset <- sce[,grepl("LAMP5", sce$annotation)]

unique(sce.subset$celltype)
# [1] "Inhib"

unique(sce.subset$annotation)
# [1] "Inh_LAMP5_2" "Inh_LAMP5_1"

# find marker genes
marker.info <- marker.info <- findMarkers(sce.subset, sce.subset$annotation, test="t", direction="up")
marker.info

# get cluster id for loop
cluster_id <- unique(levels(as.factor(sce.subset$annotation)))

# loop through cluster ids
for (i in 1:length(cluster_id)) {
    
    # just look at top 5 genes for one spatial cluster as a sanity check
    chosen <- marker.info[[i]]

    plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", "1v1_LAMP5", cluster_id[i])
    
    #make sure folder exists, if not, make it
    if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
    }
    
    # output for keeping track with loop
    print(paste("Now printing marker genes for:", cluster_id[i]))
    
    
    # visualize mean AUC
    top10 <- c(rownames(chosen[1:10,]))
    p <- scater::plotExpression(sce.subset, features=top10, x='annotation', colour_by='annotation', show_violin=TRUE, show_median=TRUE, ncol=2)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_Top10.pdf"))
    print(p)
    dev.off()
    
}


# === Volano plots ===

marker.info <- findMarkers(sce.subset, sce.subset$annotation, test="t")

# Loop through each unique cluster ID
for(cluster_id in unique(sce.subset$annotation)){

    # set plot dir
    plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", "1v1_LAMP5", cluster_id)

    # Filter the all_markers data frame for the current cluster
    cluster_markers <- as.data.frame(marker.info[[cluster_id]])

    volcano_plot <- EnhancedVolcano(
        cluster_markers,
        lab = row.names(cluster_markers),
        x = 'summary.logFC',
        y = 'p.value',
        title = paste('Volcano plot of Cluster', cluster_id, 'markers'),
        FCcutoff = 2,
        pCutoff = 10e-2,
        drawConnectors = TRUE,
        widthConnectors = 0.75
    )
    
    # Optionally, save the volcano plot to a file
    ggsave(filename = here(plot_dir, paste0("volcano_plot_cluster_", cluster_id, ".png")), plot = volcano_plot)
    
}



# ============= VIP 1 vs VIP 2 =============

# run find markers just comparing VIP 1 and 2
sce.subset <- sce[,grepl("VIP", sce$annotation)]

unique(sce.subset$celltype)
# [1] "Inhib"

unique(sce.subset$annotation)
# [1] "Inh_Vip_1" "Inh_Vip_2"

# find marker genes
marker.info <- marker.info <- findMarkers(sce.subset, sce.subset$annotation, test="t", direction="up")
marker.info

# get cluster id for loop
cluster_id <- unique(levels(as.factor(sce.subset$annotation)))

# loop through cluster ids
for (i in 1:length(cluster_id)) {
    
    # just look at top 5 genes for one spatial cluster as a sanity check
    chosen <- marker.info[[i]]

    plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", "1v1_VIP", cluster_id[i])
    
    #make sure folder exists, if not, make it
    if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
    }
    
    # output for keeping track with loop
    print(paste("Now printing marker genes for:", cluster_id[i]))
    
    
    # visualize mean AUC
    top10 <- c(rownames(chosen[1:10,]))
    p <- scater::plotExpression(sce.subset, features=top10, x='annotation', colour_by='annotation', show_violin=TRUE, show_median=TRUE, ncol=2)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_Top10.pdf"))
    print(p)
    dev.off()
    
}


# === Volano plots ===

marker.info <- findMarkers(sce.subset, sce.subset$annotation, test="t")

# Loop through each unique cluster ID
for(cluster_id in unique(sce.subset$annotation)){

    # set plot dir
    plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", "1v1_VIP", cluster_id)

    # Filter the all_markers data frame for the current cluster
    cluster_markers <- as.data.frame(marker.info[[cluster_id]])

    volcano_plot <- EnhancedVolcano(
        cluster_markers,
        lab = row.names(cluster_markers),
        x = 'summary.logFC',
        y = 'p.value',
        title = paste('Volcano plot of Cluster', cluster_id, 'markers'),
        FCcutoff = 2,
        pCutoff = 10e-2,
        drawConnectors = TRUE,
        widthConnectors = 0.75
    )
    
    # Optionally, save the volcano plot to a file
    ggsave(filename = here(plot_dir, paste0("volcano_plot_cluster_", cluster_id, ".png")), plot = volcano_plot)
    
}

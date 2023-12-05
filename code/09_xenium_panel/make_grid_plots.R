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
# [1] "sample_id"              "in_tissue"              "array_row"             
# [4] "array_col"              "10x_graphclust"         "10x_kmeans_10_clusters"
# [7] "10x_kmeans_2_clusters"  "10x_kmeans_3_clusters"  "10x_kmeans_4_clusters" 
# [10] "10x_kmeans_5_clusters"  "10x_kmeans_6_clusters"  "10x_kmeans_7_clusters" 
# [13] "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                   
# [16] "sum_umi"                "sum_gene"               "expr_chrM"             
# [19] "expr_chrM_ratio"        "ManualAnnotation"       "slide"                 
# [22] "array"                  "brnum"                  "species"               
# [25] "replicate"              "overlaps_tissue"        "sum"                   
# [28] "detected"               "subsets_mito_sum"       "subsets_mito_detected" 
# [31] "subsets_mito_percent"   "total"                  "qc_lib_size"           
# [34] "qc_mito"                "qc_detected"            "discard"               
# [37] "sizeFactor"             "row"                    "col"                   
# [40] "cluster.init"           "spatial.cluster"     


# Let's start with broad marker genes

broad_markers <- c("SNAP25", "SYT1",     # neurons
             'SLC17A7', "SLC17A6", # excitatory neurons
             "GAD1", "GAD2",       # inhibitory neurons
             "MBP", "MOBP",        # oligodendrocytes
             "GFAP", "AQP4",       # astrocytes
             "C3", "CSF1R",        # microglia
             "PDGFRA", "VCAN",     # OPCs
             "CLDN5", "FLT1"      # endothelial cells
            )
             
inhibitory_markers <- c("SST", "CORT", "CRHBP", "NPY", "CHODL", "NOS1", "TAC1", "TACR1",
                        "PVALB", 
                        "VIP",
                        "LAMP5",
                        "PRKCD", "PENK", "PDYN", "CRH", "NTS","CCK"
                        )
             

amy_markers <- c("PRKCB", "CYP26B1", "HPCAL1","CAPS", "RGS4", "NEFM")


# ============ plot grids ============

# plot broad celltype markers
for (i in 1:length(broad_markers)){
    p_list <- vis_grid_gene(
        spe,
        geneid= broad_markers[i],
        spatial = FALSE,
        auto_crop = TRUE,
        return_plots = TRUE,
        pdf_file = NULL,
        )
    plot_list_reordered <- c(p_list[5], p_list[6], p_list[4],
                             p_list[7], p_list[8], p_list[3],
                             p_list[1], p_list[2])
    cowplot::plot_grid(plotlist = plot_list_reordered, ncol = 3)
    ggsave(here(plot_dir, "broad_markers", paste0(broad_markers[i], ".pdf")), width = 20, height = 20)
}

# plot inhibitory celltype markers
for (i in 1:length(inhibitory_markers)){
    p_list <- vis_grid_gene(
        spe,
        geneid= inhibitory_markers[i],
        spatial = FALSE,
        auto_crop = TRUE,
        return_plots = TRUE,
        pdf_file = NULL,
    )
    plot_list_reordered <- c(p_list[5], p_list[6], p_list[4],
                             p_list[7], p_list[8], p_list[3],
                             p_list[1], p_list[2])
    cowplot::plot_grid(plotlist = plot_list_reordered, ncol = 3)
    ggsave(here(plot_dir, "inhibitory_types", paste0(inhibitory_markers[i], ".pdf")), width = 20, height = 20)
}

# plot amygdala subregion markers
for (i in 1:length(amy_markers)){
    p_list <- vis_grid_gene(
        spe,
        geneid= amy_markers[i],
        spatial = FALSE,
        auto_crop = TRUE,
        return_plots = TRUE,
        pdf_file = NULL,
    )
    plot_list_reordered <- c(p_list[5], p_list[6], p_list[4],
                             p_list[7], p_list[8], p_list[3],
                             p_list[1], p_list[2])
    cowplot::plot_grid(plotlist = plot_list_reordered, ncol = 3)
    ggsave(here(plot_dir, "amygdala_subregions", paste0(amy_markers[i], ".pdf")), width = 20, height = 20)
}

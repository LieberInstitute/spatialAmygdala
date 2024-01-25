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


load(here("processed-data", "02_build_spe", "spe_raw.Rdata"), verbose = TRUE)
spe.2 <- spe



# Drop any duplicate coolData
colData(spe) <- colData(spe)[ , !duplicated(colnames(colData(spe)))]
colnames(colData(spe))



# ==== combine SPEs ====

# get common genes
common_genes <- intersect(rownames(spe.1), rownames(spe.2))
spe.1 <- spe.1[common_genes,]
spe.2 <- spe.2[common_genes,]

# combine
spe <- cbind(spe.1, spe.2)

# drop out of tissue spots
spe <- spe[, spe$in_tissue]

# visualize some QC metrics
unique(colnames(colData(spe)))

unique(spe$brnum)
# [1] Br9469      Br6471      85v_AMY_SVB

# replace 85v_AMY_SVB with Br6471
spe$brnum[spe$brnum == "85v_AMY_SVB"] <- "Br6471"


# normalize
library(scuttle)
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)




# Let's start with broad marker genes

# read in xlsx file from processed_dir
custom_markers <- read.csv(here(processed_dir, "Amygdala_Xenium_panel_2.0_final.csv"), header = TRUE, stringsAsFactors = FALSE)


# read csv
SVGs <- read.csv(here("processed-data","05_feature_selection", "nnSVG_summary.csv"), header = TRUE, stringsAsFactors = FALSE)
#   gene_id gene_name gene_type overall_rank average_rank n_withinTop100
# 1     MBP       MBP      gene          1.5        2.000              8
# 2     AVP       AVP      gene          1.5        2.000              1
# 3    NGFR      NGFR      gene          3.0        5.000              1
# 4    ENC1      ENC1      gene          4.0        5.750              8
# 5  SNAP25    SNAP25      gene          5.0        7.625              8
# 6    GFAP      GFAP      gene          6.0        7.875              8


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
# custom markers
for (i in 1:length(custom_markers$Gene)){
    p_list <- vis_grid_gene(
        spe,
        geneid= custom_markers$Gene[i],
        spatial = FALSE,
        auto_crop = TRUE,
        return_plots = TRUE,
        pdf_file = NULL,
    )
    plot_list_reordered <- c(p_list[5], p_list[6], p_list[4],
                             p_list[7], p_list[8], p_list[3],
                             p_list[1], p_list[2])
    cowplot::plot_grid(plotlist = plot_list_reordered, ncol = 3)
    ggsave(here(plot_dir, "custom_markers", paste0(custom_markers$Region[i],"_", custom_markers$Gene[i], ".pdf")), width = 20, height = 20)
}

# nnSVGs
for (i in 1:length(SVGs$gene_id)){
    p_list <- vis_grid_gene(
        spe,
        geneid= SVGs$gene_id[i],
        spatial = FALSE,
        auto_crop = TRUE,
        return_plots = TRUE,
        pdf_file = NULL,
    )
    plot_list_reordered <- c(p_list[5], p_list[6], p_list[4],
                             p_list[7], p_list[8], p_list[3],
                             p_list[1], p_list[2])
    cowplot::plot_grid(plotlist = plot_list_reordered, ncol = 3)
    ggsave(here(plot_dir, "nnSVGs", paste0("Rank", SVGs$overall_rank[i],"_", SVGs$gene_id[i], ".pdf")), width = 20, height = 20)
}

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




# ======== plotting all brains ======

# ========= Spotplots ============
#

rm(spe.1)
rm(spe.2)
#spe.amy <-spe
spe <- spe.amy

# drop unused levels
spe$brnum <- droplevels(spe$brnum)

for (i in 1:length(unique(spe.amy$brnum))) {
    
    print(i)
    #subset to brnum only
    brain <- as.character(unique(spe.amy$brnum)[[i]])
    spe <- spe.amy[,spe.amy$brnum == brain]
    
    for (j in 1:length(custom_markers$Gene)){
        p_list <- vis_grid_gene(
            spe,
            geneid= custom_markers$Gene[j],
            spatial = FALSE,
            auto_crop = TRUE,
            return_plots = TRUE,
            pdf_file = NULL,
        )
        plot_list_reordered <- c(p_list[5], p_list[6], p_list[4],
                                 p_list[7], p_list[8], p_list[3],
                                 p_list[1], p_list[2])
        cowplot::plot_grid(plotlist = plot_list_reordered, ncol = 3)
        ggsave(here(plot_dir, "custom_markers", brain,paste0(custom_markers$Region[j],"_", custom_markers$Gene[j], ".pdf")), width = 20, height = 20)
    }
}
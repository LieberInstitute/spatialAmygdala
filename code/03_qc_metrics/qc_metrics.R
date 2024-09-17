library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(scater)
library(patchwork)
library(ggpubr)
library(SpotSweeper)


# Save directories
plot_dir = here("plots", "03_qc_metrics")
processed_dir = here("processed-data","03_qc_metrics")


# ===== First round of samples (#1) =====
load(here("processed-data", "02_build_spe", "spe_raw-1st.Rdata"), verbose = TRUE)
spe.1 <- spe
spe.1
# class: SpatialExperiment 
# dim: 28412 39936 
# metadata(0):
#     assays(1): counts
# rownames(28412): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(6): source type ... gene_name Symbol.uniq
# colnames(39936): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ...
# TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
# colData names(26): sample_id in_tissue ... replicate overlaps_tissue
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
#     spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

unique(spe.1$brnum)
# [1] Br8325
# Levels: Br8325

unique(spe.1$sample_id)
# [1] "V13M06-387_A1" "V13M06-387_B1" "V13M06-387_C1" "V13M06-387_D1"
# [5] "V13M06-388_A1" "V13M06-388_B1" "V13M06-388_C1" "V13M06-388_D1"


# ===== Second round of samples (#2-3) =====
load(here("processed-data", "02_build_spe", "spe_raw-2nd.Rdata"), verbose = TRUE)
spe.2 <- spe
spe.2
# class: SpatialExperiment 
# dim: 28906 79872 
# metadata(0):
#     assays(1): counts
# rownames(28906): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(6): source type ... gene_name Symbol.uniq
# colnames(79872): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ...
# TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
# colData names(26): sample_id in_tissue ... replicate overlaps_tissue
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
#     spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

unique(spe.2$brnum)
# [1] Br9469      Br6471      85v_AMY_SVB
# 11 Levels: 65v_AMY_SVB 66v_AMY_SVB 67v_AMY_SVB 68v_AMY_SVB ... Br9469

unique(spe.2$sample_id)
# [1] "V13F27-349_A1" "V13F27-349_B1" "V13F27-349_C1" "V13F27-349_D1"
# [5] "V13F27-354_A1" "V13F27-354_B1" "V13F27-354_C1" "V13F27-354_D1"
# [9] "V13F27-359_A1" "V13F27-359_B1" "V13F27-359_C1" "V13F27-359_D1"
# [13] "V13F27-366_A1" "V13F27-366_B1" "V13F27-366_C1" "V13F27-366_D1"

# ===== Third round of samples (#4-5) =====
load(here("processed-data", "02_build_spe", "spe_raw-3rd.Rdata"), verbose = TRUE)
spe.3 <- spe
spe.3
# class: SpatialExperiment 
# dim: 29335 74880 
# metadata(0):
#     assays(1): counts
# rownames(29335): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(6): source type ... gene_name Symbol.uniq
# colnames(74880): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ...
# TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
# colData names(26): sample_id in_tissue ... replicate overlaps_tissue
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
#     spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

unique(spe.3$brnum)
# [1] Br2743 Br6423
# 11 Levels: 65v_AMY_SVB 66v_AMY_SVB 67v_AMY_SVB 68v_AMY_SVB ... Br9469

# ===== Fourth round of samples (#6) =====
load(here("processed-data", "02_build_spe", "spe_raw-4th.Rdata"), verbose = TRUE)
spe.4 <- spe
spe.4

unique(spe.4$brnum)
#[1] Br6660



# ===== Fifth round of samples (#6) =====
load(here("processed-data", "02_build_spe", "spe_raw-5th.Rdata"), verbose = TRUE)
spe.5 <- spe
spe.5
# class: SpatialExperiment 
# dim: 28637 109824 
# metadata(0):
# assays(1): counts
# rownames(28637): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(6): source type ... gene_name Symbol.uniq
# colnames(109824): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ...
#   TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
# colData names(26): sample_id in_tissue ... replicate overlaps_tissue
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

unique(spe.5$brnum)
# [1] Br9206 Br9017 Br9280 Br9192 BR9192
# 16 Levels: 65v_AMY_SVB 66v_AMY_SVB 67v_AMY_SVB 68v_AMY_SVB ... Br9469

# =========== Merging SCE objects ===========
# get common genes
common_genes <- Reduce(intersect, list(rownames(spe.1), rownames(spe.2), rownames(spe.3), rownames(spe.4), rownames(spe.4), rownames(spe.5)))

spe.1 <- spe.1[common_genes,]
spe.2 <- spe.2[common_genes,]
spe.3 <- spe.3[common_genes,]
spe.4 <- spe.4[common_genes,]
spe.5 <- spe.5[common_genes,]

# combine
spe <- Reduce(cbind, list(spe.1, spe.2, spe.3,spe.4, spe.5))

# drop out of tissue spots
spe <- spe[, spe$in_tissue]

# visualize some QC metrics
unique(colnames(colData(spe)))
# [1] "sample_id"              "in_tissue"              "array_row"             
# [4] "array_col"              "10x_graphclust"         "10x_kmeans_10_clusters"
# [7] "10x_kmeans_2_clusters"  "10x_kmeans_3_clusters"  "10x_kmeans_4_clusters" 
# [10] "10x_kmeans_5_clusters"  "10x_kmeans_6_clusters"  "10x_kmeans_7_clusters" 
# [13] "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                   
# [16] "sum_umi"                "sum_gene"               "expr_chrM"             
# [19] "expr_chrM_ratio"        "ManualAnnotation"       "slide"                 
# [22] "array"                  "brnum"                  "species"               
# [25] "replicate"              "overlaps_tissue"  


unique(spe$brnum)
# [1] Br8325      Br9469      Br6471      85v_AMY_SVB Br2743      Br6423     
# 11 Levels: Br8325 65v_AMY_SVB 66v_AMY_SVB 67v_AMY_SVB ... Br9469

# replace 85v_AMY_SVB with Br6471
spe$brnum[spe$brnum == "85v_AMY_SVB"] <- "Br6471"
spe$brnum[spe$brnum == "Br9192"] <- "BR9192"


# reset levels
spe$brnum <- factor(spe$brnum)

# check
unique(spe$brnum)
# [1] Br8325 Br9469 Br6471 Br2743 Br6423 Br6660 Br9206 Br9017 Br9280 BR9192
# 10 Levels: Br8325 Br2743 Br6423 Br6471 Br6660 Br9469 Br9017 BR9192 ... Br9280

length(unique(spe$sample_id))
# [1] 69

spe
# class: SpatialExperiment 
# dim: 25814 309099 
# metadata(0):
# assays(1): counts
# rownames(25814): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(6): source type ... gene_name Symbol.uniq
# colnames(309099): AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 ...
#   TTGTTTCCATACAACT-1 TTGTTTGTGTAAATTC-1
# colData names(26): sample_id in_tissue ... replicate overlaps_tissue
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

# save combined object as rds
saveRDS(spe, here(processed_dir, "spe_combined_noQC.Rds"))


#  ========= Violin plots of qc metrics across samples ==========

# mito ratio
png(width=15, height=5, here(plot_dir,"Violin_mito_ratio.png"), res=300, units="in")
p1 <- plotColData(spe, x="sample_id", y="expr_chrM_ratio", colour_by="brnum") + 
    #scale_y_log10() + 
    ggtitle("Mitochondrial Percent") +
    geom_hline(aes(yintercept = 0.3), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$brnum, 
               scales = "free_x",
               switch = "x",
               nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")

p1
dev.off()

# sum umi
png(width=15, height=5, here(plot_dir,"Violin_sum_umi.png"), res=300, units="in")
p2 <- plotColData(spe, x="sample_id", y="sum_umi", colour_by="brnum") + 
    scale_y_log10() + 
    ggtitle("Sum UMI") +
    geom_hline(aes(yintercept = 1000), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$brnum, 
               scales = "free_x",
               switch = "x",
               nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")

p2
dev.off()

# sum gene
png(width=15, height=5, here(plot_dir,"Violin_sum_gene.png"), res=300, units="in")
p3 <- plotColData(spe, x="sample_id", y="sum_gene", colour_by="brnum") + 
    scale_y_log10() + 
    ggtitle("Sum Gene") +
    geom_hline(aes(yintercept = 1000), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$brnum, 
               scales = "free_x",
               switch = "x",
               nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")

p3
dev.off()

png(width=15, height=10, here(plot_dir,"Violin_all_metrics.png"), res=300, units="in")
p1 <- p1 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank()
                 )

p2 <- p2 + theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()
                )

p1/p2/p3
dev.off()





 # ========= Spotplots ============
#
spe
spe.amy <- spe
unique(spe.amy$brnum)

# drop unused levels
spe.amy$brnum <- droplevels(spe.amy$brnum)

for (i in 1:length(unique(spe.amy$brnum))) {

    print(i)
    #subset to brnum only
    brain <- as.character(unique(spe.amy$brnum)[[i]])
    spe <- spe.amy[,spe.amy$brnum == brain]
    
    # mito percent
    png(width=20, height=10, here(plot_dir,"Spotplots","mito_percent", paste0("MitoPerc_", brain,".png")), res=300, units="in")
    sample_ids <- unique(spe$sample_id)
    plots <- vis_grid_gene(
        spe = spe,
        geneid = "expr_chrM_ratio",
        point_size=1.5,
        return_plots=TRUE,
        assay="counts"
      )
    print(cowplot::plot_grid(plotlist = plots, ncol = 4))
    dev.off()
    
    # libary size
    png(width=20, height=10, here(plot_dir,"Spotplots","library_size", paste0("Umi_", brain,".png")), res=300, units="in")
    sample_ids <- unique(spe$sample_id)
    plots <- vis_grid_gene(
        spe = spe,
        geneid = "sum_umi",
        point_size=1.5,
        return_plots=TRUE,
        assay="counts"
      )
    print(cowplot::plot_grid(plotlist = plots, ncol = 4))
    dev.off()
    
    # unique genes
    png(width=20, height=10, here(plot_dir,"Spotplots","unique_genes", paste0("Genes_", brain,".png")), res=300, units="in")
    sample_ids <- unique(spe$sample_id)
    plots <- vis_grid_gene(
        spe = spe,
        geneid = "sum_gene",
        point_size=1.5,
        return_plots=TRUE,
        assay="counts"
      )
    print(cowplot::plot_grid(plotlist = plots, ncol = 4))
    dev.off()

}

# =============== Calculate QC Metrics =================
#
#
# will wait to use SpotSweeper rather than waste time on hard thresholds
#
#
#
#

library(SpotSweeper)

spe <- spe.amy
rm(spe.amy)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
#FALSE  TRUE 
#28399    13 

rowData(spe)$gene_name[is_mito]
# [1] "MT-ND1"  "MT-ND2"  "MT-CO1"  "MT-CO2"  "MT-ATP8" "MT-ATP6" "MT-CO3"  "MT-ND3"  "MT-ND4L" "MT-ND4"  "MT-ND5"  "MT-ND6"  "MT-CYB" 


# ======= SpotSweeper ======
# library size
spe <- localOutliers(spe, metric="sum_umi",direction="lower", log=TRUE)

# unique genes
spe <- localOutliers(spe, metric="sum_gene", direction="lower", log=TRUE)

# mitochondrial percent
spe <- localOutliers(spe, metric="expr_chrM_ratio", direction="higher", log=FALSE)


# combine all outliers into "local_outliers" column
spe$local_outliers <- as.logical(spe$sum_umi_outliers) | 
    as.logical(spe$sum_gene_outliers) | 
    as.logical(spe$expr_chrM_ratio_outliers)

# save spe
save(spe, file = here(processed_dir, "spe_local_outliers.Rdata"))

# load
load(here(processed_dir, "spe_local_outliers.Rdata"))
spe


# ======= Visualzing QC =========

table(spe$local_outliers)
# FALSE   TRUE 
#207745    867 

# percent
table(spe$local_outliers)/ncol(spe)*100
#      FALSE       TRUE 
# 99.5843959  0.4156041 

# subset by brnum, then plot pdf of all samples
for (i in 1:length(unique(spe$brnum))) {
    brain <- as.character(unique(spe$brnum)[[i]])
    spe.tmp <- spe[,spe$brnum == brain]
    

    plotQCpdf(spe.tmp,
                    metric="sum_umi_log",
                    outliers="local_outliers",
                    point_size=2,
                    stroke=0.8,
                    fname=here(plot_dir,"Spotplots", "local_outliers", paste0("LocalOutliers_", brain, ".pdf")))

}

# number of local outliers per brnum
table(spe$brnum, spe$local_outliers)
#          FALSE  TRUE
#   Br8325 29776   109
#   Br2743 32691   164
#   Br6423 34222   131
#   Br6471 37602   193
#   Br6660 36884   201
#   Br9469 36570    69


#  ======= Violin plots ========

# mito ratio
png(width=15, height=5, here(plot_dir,"Violin_mito_ratio_z.png"), res=300, units="in")
p1 <- plotColData(spe, x="sample_id", y="expr_chrM_ratio_z", colour_by="local_outliers") + 
    ggtitle("Mitochondrial Percent") +
    geom_hline(aes(yintercept = 3), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$brnum, 
               scales = "free_x",
               switch = "x",
               nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")

p1
dev.off()

# sum umi
png(width=15, height=5, here(plot_dir,"Violin_sum_umi_z.png"), res=300, units="in")
p2 <- plotColData(spe, x="sample_id", y="sum_umi_z", colour_by="local_outliers") + 
    ggtitle("Sum UMI") +
    geom_hline(aes(yintercept = -3), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$brnum, 
               scales = "free_x",
               switch = "x",
               nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")

p2
dev.off()

# sum gene
png(width=15, height=5, here(plot_dir,"Violin_sum_gene_z.png"), res=300, units="in")
p3 <- plotColData(spe, x="sample_id", y="sum_gene_z", colour_by="local_outliers") + 
    ggtitle("Sum Gene") +
    geom_hline(aes(yintercept = -3), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$brnum, 
               scales = "free_x",
               switch = "x",
               nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")

p3
dev.off()

png(width=15, height=10, here(plot_dir,"Violin_all_metrics_z.png"), res=300, units="in")
p1 <- p1 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank()
)

p2 <- p2 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank()
)

p1/p2/p3
dev.off()
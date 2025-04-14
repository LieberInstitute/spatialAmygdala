library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(scater)
library(patchwork)
library(ggpubr)
library(SpotSweeper)


# Save directories
plot_dir = here("plots", "Visium", "03_qc_metrics", "stitched")
processed_dir = here("processed-data","Visium","03_qc_metrics")

# save combined object as rds
spe <- readRDS(here(processed_dir, "spe_stitched_combined_noQC.Rds"))
spe

# ========= Add QC metrics to spe ==========
colnames(colData(spe))
#  [1] "sample_id"                   "in_tissue"                  
#  [3] "X10x_graphclust"             "X10x_kmeans_10_clusters"    
#  [5] "X10x_kmeans_2_clusters"      "X10x_kmeans_3_clusters"     
#  [7] "X10x_kmeans_4_clusters"      "X10x_kmeans_5_clusters"     
#  [9] "X10x_kmeans_6_clusters"      "X10x_kmeans_7_clusters"     
# [11] "X10x_kmeans_8_clusters"      "X10x_kmeans_9_clusters"     
# [13] "key"                         "sum_umi"                    
# [15] "sum_gene"                    "expr_chrM"                  
# [17] "expr_chrM_ratio"             "ManualAnnotation"           
# [19] "capture_area"                "group"                      
# [21] "barcode"                     "array_row_original"         
# [23] "array_col_original"          "array_row"                  
# [25] "array_col"                   "pxl_col_in_fullres_rounded" 
# [27] "pxl_row_in_fullres_rounded"  "pxl_row_in_fullres_original"
# [29] "pxl_col_in_fullres_original" "overlap_key"                
# [31] "exclude_overlapping"  

# QC metrics are already there, so we won't add them again

#  ========= Violin plots of qc metrics across samples ==========
# mito ratio
png(width=15, height=5, here(plot_dir,"Violin_mito_ratio.png"), res=300, units="in")
p1 <- plotColData(spe, x="capture_area", y="expr_chrM_ratio", colour_by="sample_id") + 
    #scale_y_log10() + 
    ggtitle("Mitochondrial Percent") +
    geom_hline(aes(yintercept = 0.3), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$sample_id, 
                scales = "free_x",
                switch = "x",
                nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")

p1
dev.off()

# sum umi
png(width=15, height=5, here(plot_dir,"Violin_sum_umi.png"), res=300, units="in")
p2 <- plotColData(spe, x="capture_area", y="sum_umi", colour_by="sample_id") + 
    scale_y_log10() + 
    ggtitle("Sum UMI") +
    geom_hline(aes(yintercept = 1000), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$sample_id, 
               scales = "free_x",
               switch = "x",
               nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")

p2
dev.off()

# sum gene
png(width=15, height=5, here(plot_dir,"Violin_sum_gene.png"), res=300, units="in")
p3 <- plotColData(spe, x="capture_area", y="sum_gene", colour_by="sample_id") + 
    scale_y_log10() + 
    ggtitle("Sum Gene") +
    geom_hline(aes(yintercept = 1000), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$sample_id, 
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
spe <- spe[,spe$in_tissue]

plotQCpdf(spe,
        metric="sum_umi",
        outliers=NULL,
        point_size=.7,
        stroke=0.8,
        fname=here(plot_dir,"Spotplots", paste0("sum_umi_Spotplots.pdf")))

plotQCpdf(spe,
        metric="sum_gene",
        outliers=NULL,
        point_size=.7,
        stroke=0.8,
        fname=here(plot_dir,"Spotplots", paste0("sum_gene_Spotplots.pdf")))

plotQCpdf(spe,
        metric="expr_chrM_ratio",
        outliers=NULL,
        point_size=.7,
        stroke=0.8,
        fname=here(plot_dir,"Spotplots", paste0("Mito_Spotplots.pdf")))


# =============== Calculate QC Metrics =================

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
save(spe, file = here(processed_dir, "spe_stitched_local_outliers.Rdata"))



# load
load(here(processed_dir, "spe_stitched_local_outliers.Rdata"))
spe


# ======= Visualzing QC =========

plotQCpdf(spe,
        metric="sum_umi",
        outliers="local_outliers",
        point_size=.7,
        stroke=0.8,
        fname=here(plot_dir,"Spotplots", paste0("sum_umi_Spotplots_outliers.pdf")))

plotQCpdf(spe,
        metric="sum_gene",
        outliers="local_outliers",
        point_size=.7,
        stroke=0.8,
        fname=here(plot_dir,"Spotplots", paste0("sum_gene_Spotplots_outliers.pdf")))

plotQCpdf(spe,
        metric="expr_chrM_ratio",
        outliers="local_outliers",
        point_size=.7,
        stroke=0.8,
        fname=here(plot_dir,"Spotplots", paste0("Mito_Spotplots_outliers.pdf")))



# number of local outliers per sample_id
table(spe$sample_id, spe$local_outliers)
#          FALSE  TRUE
#   Br2743 32665   190
#   Br6423 34217   136
#   Br6471 37578   217
#   Br6660 36860   225
#   Br8325 29753   132
#   Br9017 26335    47
#   Br9192 28553    55
#   Br9206 26798    52
#   Br9280 26585   137
#   Br9469 36577    62

#  ======= Violin plots ========

# mito ratio
png(width=15, height=5, here(plot_dir,"Violin_expr_chrM_ratio_z.png"), res=300, units="in")
p1 <- plotColData(spe, x="capture_area", y="expr_chrM_ratio_z", colour_by="local_outliers") + 
    ggtitle("Mitochondrial Percent") +
    geom_hline(aes(yintercept = 3), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$sample_id, 
               scales = "free_x",
               switch = "x",
               nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")

p1
dev.off()

# sum umi
png(width=15, height=5, here(plot_dir,"Violin_sum_umi_z.png"), res=300, units="in")
p2 <- plotColData(spe, x="capture_area", y="sum_umi_z", colour_by="local_outliers") + 
    ggtitle("Sum UMI") +
    geom_hline(aes(yintercept = -3), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$sample_id, 
               scales = "free_x",
               switch = "x",
               nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")

p2
dev.off()

# sum gene
png(width=15, height=5, here(plot_dir,"Violin_sum_gene_z.png"), res=300, units="in")
p3 <- plotColData(spe, x="capture_area", y="sum_gene_z", colour_by="local_outliers") + 
    ggtitle("Sum Gene") +
    geom_hline(aes(yintercept = -3), linetype="dashed", color = "red") +
    #coord_flip() + 
    facet_wrap(~spe$sample_id, 
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
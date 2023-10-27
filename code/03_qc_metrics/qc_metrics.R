## qrsh -l mem_free=30G,h_vmem=30G -now n

library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(scater)
library(patchwork)


# Save directories
plot_dir = here("plots", "03_qc_metrics")
processed_dir = here("processed-data","03_qc_metrics")


# Load AMY data
load(here("processed-data", "02_build_spe", "spe_raw.Rdata"), verbose = TRUE)
spe
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

# drop out of tissue spots
spe <- spe[, spe$in_tissue]

# visualize some QC metrics

# expr_chrM
pdf(width=20, height=10, here(plot_dir,"Spotplot_expr_chrM.pdf"))
sample_ids <- unique(spe$sample_id)
plots <- lapply(sample_ids[1:8], function(sampleid) {
  vis_gene(
    spe = spe,
    geneid = "expr_chrM",
    sampleid = sampleid
  )
})
(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]] | plots[[7]] | plots[[8]])
dev.off()

# umi
pdf(width=20, height=10, here(plot_dir,"Spotplot_sum_umi.pdf"))
sample_ids <- unique(spe$sample_id)
plots <- lapply(sample_ids[1:8], function(sampleid) {
  vis_gene(
    spe = spe,
    geneid = "sum_umi",
    sampleid = sampleid
  )
})
(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]] | plots[[7]] | plots[[8]])
dev.off()

# umi
pdf(width=20, height=10, here(plot_dir,"Spotplot_expr_chrM_ratio.pdf"))
sample_ids <- unique(spe$sample_id)
plots <- lapply(sample_ids[1:8], function(sampleid) {
  vis_gene(
    spe = spe,
    geneid = "expr_chrM_ratio",
    sampleid = sampleid
  )
})
(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]] | plots[[7]] | plots[[8]])
dev.off()

# umi
pdf(width=20, height=10, here(plot_dir,"Spotplot_sum_gene.pdf"))
sample_ids <- unique(spe$sample_id)
plots <- lapply(sample_ids[1:8], function(sampleid) {
  vis_gene(
    spe = spe,
    geneid = "sum_gene",
    sampleid = sampleid
  )
})
(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]] | plots[[7]] | plots[[8]])
dev.off()


# =============== Calculate QC Metrics =================
# to perform standard QC in similar fashion to snRNA-seq data we'll use the scater package

library(scater)
library(ggspavis)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
#FALSE  TRUE 
#28399    13 

rowData(spe)$gene_name[is_mito]
# [1] "MT-ND1"  "MT-ND2"  "MT-CO1"  "MT-CO2"  "MT-ATP8" "MT-ATP6" "MT-CO3"  "MT-ND3"  "MT-ND4L" "MT-ND4"  "MT-ND5"  "MT-ND6"  "MT-CYB" 


# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
colnames(colData(spe))
# [1] "sample_id"              "in_tissue"              "array_row"              "array_col"              "10x_graphclust"         "10x_kmeans_10_clusters"
# [7] "10x_kmeans_2_clusters"  "10x_kmeans_3_clusters"  "10x_kmeans_4_clusters"  "10x_kmeans_5_clusters"  "10x_kmeans_6_clusters"  "10x_kmeans_7_clusters" 
# [13] "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                    "sum_umi"                "sum_gene"               "expr_chrM"             
# [19] "expr_chrM_ratio"        "ManualAnnotation"       "subject"                "region"                 "sex"                    "age"                   
# [25] "diagnosis"              "sample_id_complete"     "count"                  "sum"                    "detected"               "subsets_mito_sum"      
# [31] "subsets_mito_detected"  "subsets_mito_percent"   "total"    


# ====================== Library Size ========================

# histogram of library sizes
pdf(here(plot_dir, "Histogram_LibrarySize.pdf"))
hist(colData(spe)$sum, breaks = 80, xlim=c(0,15000))
dev.off()

# histogram of mitochondrial read proportions
pdf(here(plot_dir, "Histogram_MitoPercent.pdf"))
hist(colData(spe)$subsets_mito_percent, breaks = 20)
dev.off()


# get library size cutoff
qc_lib_size <- colData(spe)$sum < 50 # this is arbitrary
table(qc_lib_size)

# add library size info
colData(spe)$qc_lib_size <- qc_lib_size

# check spatial pattern of discarded spots
pdf(here(plot_dir,"LibrarySize_vs_SumGenes.pdf"))
plotQC(spe, type = "scatter", 
       discard = "qc_lib_size",
       metric_x="sum_gene"
       )
dev.off()

# pdf(here(plot_dir,"LibrarySize_vs_Count.pdf"))
# plotQC(spe, type = "scatter", 
#        discard = "qc_lib_size",
#        metric_x="count"
# )
# dev.off()

# check spatial pattern of discarded spots
pdf(here(plot_dir,"MitoPercent_vs_SumGenes.pdf"))
plotQC(spe, type = "scatter", 
       metric_y="subsets_mito_percent",
       metric_x="sum_gene",
       threshold_x = 150,
)
dev.off()

# check spatial pattern of discarded spots
pdf(here(plot_dir,"MitoPercent_vs_SumUMI.pdf"))
plotQC(spe, type = "scatter", 
       metric_y="subsets_mito_percent",
       metric_x="sum_umi",
       threshold_x = 150,
)
dev.off()

# ====== select QC threshold for library size (total genes) ========
#FALSE  TRUE 
#28090  1795 

colData(spe)$qc_lib_size <- qc_lib_size

# Spotplots
pdf(width=20, height=10, here(plot_dir,"Spotplot_QC_lib_size_250.pdf"))
sample_ids <- unique(spe$sample_id)
plots <- lapply(sample_ids[1:8], function(sampleid) {
  vis_clus(
    spe = spe,
    clustervar = "qc_lib_size",
    sampleid = sampleid
  )
})
(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]] | plots[[7]] | plots[[8]])
dev.off()


# ====== select QC threshold for unqie detected genes ========
qc_detected <- colData(spe)$detected < 50 # this is arbitrary
table(qc_detected)
#FALSE  TRUE 
#27649  2236

colData(spe)$qc_detected <- qc_detected

# Spotplots
pdf(width=20, height=10, here(plot_dir,"Spotplot_QC_detected_100.pdf"))
sample_ids <- unique(spe$sample_id)
plots <- lapply(sample_ids[1:8], function(sampleid) {
  vis_clus(
    spe = spe,
    clustervar = "qc_detected",
    sampleid = sampleid
  )
})
(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]] | plots[[7]] | plots[[8]])
dev.off()


# ====== select QC threshold for mito percent ========
qc_mito <- colData(spe)$subsets_mito_percent > 30 # this is arbitrary
table(qc_mito)
# FALSE  TRUE 
# 29736   149 

colData(spe)$qc_mito <- qc_mito

# Spotplots
pdf(width=20, height=10, here(plot_dir,"Spotplot_QC_mito_25.pdf"))
sample_ids <- unique(spe$sample_id)
plots <- lapply(sample_ids[1:8], function(sampleid) {
  vis_clus(
    spe = spe,
    clustervar = "qc_mito",
    sampleid = sampleid
  )
})
(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]] | plots[[7]] | plots[[8]])
dev.off()

# mito percent
pdf(width=20, height=10, here(plot_dir,"Spotplot_mito_percent.pdf"))
sample_ids <- unique(spe$sample_id)
plots <- lapply(sample_ids[1:8], function(sampleid) {
  vis_gene(
    spe = spe,
    geneid = "subsets_mito_percent",
    sampleid = sampleid
  )
})
(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]] | plots[[7]] | plots[[8]])
dev.off()

# Graph based clustering
pdf(width=20, height=10, here(plot_dir,"Spotplot_10x_graphclust.pdf"))
sample_ids <- unique(spe$sample_id)
plots <- lapply(sample_ids[1:8], function(sampleid) {
  vis_clus(
    spe = spe,
    clustervar = "10x_graphclust",
    sampleid = sampleid
  )
})
(plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]] | plots[[7]] | plots[[8]])
dev.off()


# =========== Violin plots of QC metrics per sample =========

discard <- qc_lib_size | qc_detected | qc_mito
spe$discard <- discard

pdf(width=10, height=5, here(plot_dir,"Violin_QC_sum.pdf"))
plotColData(spe, x="sample_id", y="sum", colour_by="discard") + 
  scale_y_log10() + ggtitle("Total count")
dev.off()

pdf(width=10, height=5, here(plot_dir,"Violin_QC_mito.pdf"))
plotColData(spe, x="sample_id", y="subsets_mito_percent", colour_by="discard") + 
  scale_y_log10() + ggtitle("Mitochondrial Percent")
dev.off()

save(spe, file=here(processed_dir, "spe_discarded.Rdata"))

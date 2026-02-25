library(SpatialExperiment)
library(RcppML)
library(here)
library(Matrix)
library(ggplot2)
library(patchwork)
library(projectR)

plot_dir <- here("plots","98_NMF", "NMF_Yu")
processed_dir <- here("processed-data", "Visium","98_NMF")

# load NMF results
load(here(processed_dir,"NMF_Yu", "RcppML_NMF_Yu.rda"))
#x

# load Spatial object
load(here("processed-data","Visium","06_batch_correction", "spe_harmony.Rdata"))
spe
# class: SpatialExperiment 
# dim: 36601 317173 
# metadata(0):
# assays(2): counts logcounts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type gene_search
# colnames(317173): AAACAAGTATCTCCCA-1_V13Y24-346_A1
#   AAACAATCTACTAGCA-1_V13Y24-346_A1 ... TTGTTTGTATTACACG-1_V13F27-349_D1
#   TTGTTTGTGTAAATTC-1_V13F27-349_D1
# colData names(42): sample_id in_tissue ... sizeFactor slide_id
# reducedDimNames(3): PCA PCA-HARMONY_sample_slide PCA-HARMONY_sample
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

# load Single Nucleus object
load(here("processed-data", "snRNAseq", "yu_sce_gtf.rda"))
sce <- sce.amy
sce

# extract patterns
patterns <- t(x$h)
colnames(patterns) <- paste("NMF", 1:100, sep = "_")

loadings <- x$w
rownames(loadings) <- rownames(sce)



# ====== project loadings to spatial data =======
# drop any rownames in SPE not in scee
spe<- spe[rownames(spe) %in% rownames(sce),]

# drop any rownames in loadings not in spe
loadings <- loadings[rownames(loadings) %in% rownames(spe),]

# get common genes
common_genes <- intersect(rownames(spe), rownames(loadings))

# subset spe and loadings to common genes
spe <- spe[common_genes, ]
loadings <- loadings[common_genes, ]

logcounts <- logcounts(spe)
#data <- as.matrix(logcounts)

proj <- project(logcounts, loadings)
proj <- t(proj)
colnames(proj) <- paste("NMF", 1:100, sep = "_")



# add to reducedDims
reducedDim(spe, "NMF_proj") <- proj

# save
save(spe, file=here("processed-data","Visium","98_NMF", "spe_NMF_Yu.rda"))

spe.temp <- spe

# add each proj column to colData(spe)
for (i in 1:100){
    colData(spe.temp)[[paste0("NMF_",i)]] <- reducedDims(spe.temp)$NMF_proj[,i]
}

# drop duplicate colData
colData(spe.temp) <- colData(spe.temp)[ , !duplicated(colnames(colData(spe.temp)))]

for (i in 1:(dim(patterns)[2])){
    p_list <- vis_grid_gene(
        spe.temp,
        geneid= paste0("NMF_", i),
        spatial = FALSE,
        auto_crop = TRUE,
        return_plots = TRUE,
        pdf_file = NULL,
    )
    plot_list_reordered <- c(p_list[5], p_list[6], p_list[4],
                             p_list[7], p_list[8], p_list[3],
                             p_list[1], p_list[2])
    cowplot::plot_grid(plotlist = plot_list_reordered, ncol = 3)
    ggsave(here(plot_dir, "SpotPlots", paste0("NMF_", i, ".pdf")), width = 20, height = 20)
}


# ======== Heatmaps ========

# ==== Spatial domains ====

spe$spatial.cluster <- factor(spe$spatial.cluster)
levels(spe$spatial.cluster) <- c('WM.1', 'aBA', 'LA.1', 'BA', 'vmBA', 'GABA.n', 'LA.2', 'EC', 'WM.2', 'Endo')
unique(spe$spatial.cluster)

# create dataframe 
data <- data.frame(colData(spe), reducedDims(spe)$NMF_proj)

# aggregate NMF patterns across clusters. # grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("NMF", colnames(data))],
                      by=list(data$spatial.cluster), 
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]

p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column"
)

pdf(here(plot_dir, "Heatmap_NMF_spatial_domains.pdf"), width=12.5, height=5)
p1
dev.off()



# ======= Exporting top gene per factor ========

# function for getting top n genes for each pattern
top_genes <- function(W, n=10){
    top_genes <- apply(W, 2, function(x) names(sort(x, decreasing=TRUE)[1:n]))
    return(top_genes)
}

#get top genes
top15 <- top_genes(loadings, 15)
head(top15)

# add colnames
colnames(top15) <- colnames(patterns)

# export to csv
write.csv(top15, here(processed_dir, "NMF_Yu", "top15_genes_per_factor.csv"))

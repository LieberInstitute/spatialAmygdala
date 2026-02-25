library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("harmony")

# save directories
processed_dir <- here("processed-data", "Xenium", "04_dim_reduction")
plot_dir <- here("plots", "Xenium", "04_dim_reduction")

spe <- readRDS(here("processed-data","Xenium", "03_quality_control", "spe_spacetrooper_QCed.rds"))
spe
# class: SpatialExperiment 
# dim: 541 946252 
# metadata(8): polygons technology ... polygons technology
# assays(1): counts
# rownames(541): ABCC9 ADAMTS12 ... DeprecatedCodeword_0344
#   DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(946252): aaaaeomf-1 aaaajkhp-1 ... oimbboka-1 oimbcgpk-1
# colData names(31): cell_id transcript_counts ... QC_score low_qcscore
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id


# subset to only target genes probes
gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe.gex <- spe[gene_expression_idx,]
dim(spe.gex)
# [1]     366 946252


# ======== Dim reduction and Harmony batch correction ========
set.seed(1000)

# runPCA
spe.gex <- runPCA(spe.gex,  exprs_values = "cell_normcounts")

# run Harmony
spe.gex <- RunHarmony(spe.gex, group.by.vars="brnum")

# runUMAP
spe.gex <- runUMAP(spe.gex, dimred="HARMONY")

# ======= Copy reduced dims to orignal spe object =======
reducedDim(spe, "PCA") <- reducedDim(spe.gex, "PCA")
reducedDim(spe, "HARMONY") <- reducedDim(spe.gex, "HARMONY")
reducedDim(spe, "UMAP") <- reducedDim(spe.gex, "UMAP")


# ======== Visualization ========

# PCs
png(file.path(plot_dir, "Uncorrected_PCs_Brnum.png"), width = 8, height = 6, units = "in", res = 300)
plotReducedDim(spe.gex, dimred="PCA", ncomponents=4,
    colour_by="brnum")
dev.off()

png(file.path(plot_dir, "Corrected_PCs_Brnum.png"), width = 8, height = 6, units = "in", res = 300)
plotReducedDim(spe.gex, dimred="HARMONY", ncomponents=4,
    colour_by="brnum")
dev.off()

# UMAP
png(file.path(plot_dir, "Corrected_UMAP_Brnum.png"), width = 10, height = 10, units = "in", res = 300)
plotReducedDim(spe.gex, dimred="UMAP", ncomponents=2, colour_by="brnum", point_size=-.2)
dev.off()


rownames(spe.gex) <- rowData(spe.gex)$Symbol
markers <- c("MOBP","SLC17A7","GAD1","GULP1", "COL25A1","PDYN","TSHZ1", "SST")

# generate plots
plots <- list()
for (i in 1:length(markers)) {
  plots[[i]] <- plotReducedDim(spe.gex, dimred = "UMAP", colour_by = markers[i], point_size=0.6, by.assay.type="cell_normcounts") +
    scale_color_gradient(low="grey",high="red") +
    ggtitle(markers[i]) 
}

#print to png
png(here(plot_dir, "Corrected_UMAP_markers.png"), width=20, height=10, units="in", res=300)
patchwork::wrap_plots(plots, ncol=4)
dev.off()

# save
saveRDS(spe.gex, here("processed-data", "Xenium","04_dim_reduction", "sce_combined_harmonized_singlecell.rds"))
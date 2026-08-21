# modified from: https://github.com/LieberInstitute/spatialdACC/blob/main/code/11_differential_expression/nnSVG_precast_pseudobulk.R
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library("patchwork")
})


# load spe for k=9 without WM-CC
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
spe 

spe$layer <- spe$BS_k16_Semisupervised_wAI

spe_pseudo <-
    registration_pseudobulk(spe,
                            var_registration = "layer",
                            var_sample_id = "sample_id",
                            min_ncells = 10
    )

dim(spe_pseudo)

pca <- prcomp(t(assays(spe_pseudo)$logcounts))
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

# there is one WM spot that is very high in PC2 (PC2 > 300) and is an outlier
# we will remove it from the analysis

## save pseudobulked spe file
saveRDS(
    spe_pseudo,
    file = here("processed-data", "Visium", "08_marker_genes","pseudo_BS_K16_final_labels.rds")
)



# add perCellQC
spe_pseudo <- scuttle::addPerCellQC(spe_pseudo)

pdf(file = here("plots", "Visium", "08_marker_genes",
                paste0("pseudobulk_PC_BS_k16.pdf")),
    width = 10, height = 10)

print(plotPCA(
        spe_pseudo,
        colour_by = "layer",
        ncomponents = 2,
        point_size = 2,
        percentVar = metadata(spe_pseudo)$PCA_var_explained
    )
)
plotPCA(
    spe_pseudo,
    colour_by = "sample_id",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)

plotPCA(
        spe_pseudo,
        colour_by = "sum",
        ncomponents = 2,
        point_size = 2,
        percentVar = metadata(spe_pseudo)$PCA_var_explained
    )

plotPCA(
    spe_pseudo,
    colour_by = "detected",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)

vars <- getVarianceExplained(spe_pseudo,
                             variables = c("layer","sample_id", "sum", "detected")
)


plotExplanatoryVariables(vars)

dev.off()





# make supp figure
spe_pseudo$spatial_domain <- spe_pseudo$layer

p1 <- plotPCA(
    spe_pseudo,
    colour_by = "spatial_domain",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)

p2 <- plotPCA(
    spe_pseudo,
    colour_by = "sample_id",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)

p3 <- plotPCA(
    spe_pseudo,
    colour_by = "detected",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)

vars <- getVarianceExplained(spe_pseudo,
                             variables = c("spatial_domain","sample_id", "sum", "detected")
)

p4 <- plotExplanatoryVariables(vars)

png(file = here("plots",  "Visium", "08_marker_genes","pseudobulk_PC_BS_k16_wrapped.png"),
    width = 10, height = 10, unit="in", res=300)

wrap_plots(p1,p2,p3,p4,nrow=2)
dev.off()



# 6 outliers in PC1. get their sample_id, and spatial_domain, then Drop 6 highest PC1 spots. 
idx <- order(reducedDims(spe_pseudo)$PCA[,'PC01'], decreasing = TRUE)[1:6]
id <- colData(spe_pseudo)[idx, c("sample_id", "layer")]
# DataFrame with 6 rows and 2 columns
# Br9192_CHAT      Br9192     CHAT
# Br9280_HPC       Br9280     HPC 
# Br9280_AI        Br9280     AI  
# Br9280_CHAT      Br9280     CHAT
# Br2743_HPC       Br2743     HPC 
# Br2743_CHAT      Br2743     CHAT
spe_pseudo <- spe_pseudo[, -idx]

# rerun PCA
pca <- prcomp(t(assays(spe_pseudo)$logcounts))
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)


# replot


# make supp figure
spe_pseudo$spatial_domain <- spe_pseudo$layer

p1 <- plotPCA(
    spe_pseudo,
    colour_by = "spatial_domain",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)

p2 <- plotPCA(
    spe_pseudo,
    colour_by = "sample_id",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)

p3 <- plotPCA(
    spe_pseudo,
    colour_by = "detected",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)

vars <- getVarianceExplained(spe_pseudo,
                             variables = c("spatial_domain","sample_id", "sum", "detected")
)

p4 <- plotExplanatoryVariables(vars)

png(file = here("plots",  "Visium", "08_marker_genes","pseudobulk_PC_BS_k16_wrapped_dropped.png"),
    width = 10, height = 10, unit="in", res=300)

wrap_plots(p1,p2,p3,p4,nrow=2)
dev.off()


# save final pseudobulked spe file
saveRDS(
    spe_pseudo,
    file = here("processed-data", "Visium", "08_marker_genes","pseudo_BS_K16_final_labels_dropped.rds")
)   
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
    library("harmony")
    library("scater")
})

k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

processed_dir <- here("processed-data", "Xenium", "07_SEraster")
plots_dir <- here("plots", "Xenium", "07_SEraster", "BayesSpace")

# load 9017, 9192, 9206, 9280
rast_Br9017 <- readRDS(here(processed_dir, paste0("spe_rasterized_gene_expression_sum_Br9017_259um.rds")))
rast_Br9192 <- readRDS(here(processed_dir, paste0("spe_rasterized_gene_expression_sum_Br9192_259um.rds")))
rast_Br9206 <- readRDS(here(processed_dir, paste0("spe_rasterized_gene_expression_sum_Br9206_259um.rds")))
rast_Br9280 <- readRDS(here(processed_dir, paste0("spe_rasterized_gene_expression_sum_Br9280_259um.rds")))

rast_Br9017$sample_id <- "Br9017"
rast_Br9192$sample_id <- "Br9192"
rast_Br9206$sample_id <- "Br9206"
rast_Br9280$sample_id <- "Br9280"


# cbind
spe <- cbind(rast_Br9017, rast_Br9192, rast_Br9206, rast_Br9280)
spe


# ====== Normalization, PCA, Harmony ======

counts(spe) <- assay(spe, "pixelval")


set.seed(1000)
spe <- scuttle::logNormCounts(spe)
spe <- runPCA(spe, ncomponents=50, exprs_values="counts")

# plot pca
pdf(here(plots_dir, paste0("PCA_plot_55_new.pdf")), width=8, height=6)
plotPCA(spe, colour_by="sample_id") +
    ggtitle("PCA plot colored by sample_id")
dev.off()



spe <- runUMAP(spe, dimred="PCA", ncomponents=20, min_dist=0.3, name="UMAP_PCA")

pdf(here(plots_dir, paste0("UMAP_uncorrected_new.pdf")), width=8, height=6)
plotUMAP(spe, colour_by="sample_id", dimred="UMAP_PCA") +
    ggtitle("UMAP on PCA")
dev.off()



spe <- RunHarmony(spe, "sample_id")
spe <- runUMAP(spe, dimred="HARMONY", ncomponents=10, min_dist=0.3, name="UMAP-HARMONY")

# plot to check
pdf(here(plots_dir, paste0("UMAP_PCA_new.pdf")), width=12, height=5)
p1 <- plotUMAP(spe, colour_by="sample_id", dimred="UMAP_PCA") +
    ggtitle("UMAP on PCA")

p2 <- plotUMAP(spe, colour_by="sample_id", dimred="UMAP-HARMONY") +
    ggtitle("UMAP on HARMONY")
p1+p2
dev.off()



# save
saveRDS(spe, here(processed_dir, "spe_rasterized_harmony_259.rds"))
library("here")
library("SpatialExperiment")
library("SummarizedExperiment")
library("SingleCellExperiment")
library("spatialLIBD")
library("ggplot2")
library("patchwork")
library("Banksy")
library("Seurat")
library("scater")
library("scran")
library("harmony")
library("escheR")
library("data.table")


# Save directories
plot_dir = here("plots","Visium", "07_clustering","BANKSY")
processed_dir = here("processed-data","07_clustering")

spe <- readRDS(here("processed-data", "Visium","07_clustering", "BANKSY", "spe_banksy_dropped_samples.rds"))
colnames(colData(spe))


# loop through each sample_id and create a pdf on a separete page
pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(spe$clust_HARMONY_M0_lam0.4_k50_res0.4)))
pdf(file = here::here("plots", "Visium", "07_clustering", "BANKSY", "banksy_all_samples_smoothed_0.2lambda.pdf"), width = 10, height = 10)
for (id in unique(colData(spe)$sample_id)) {
    spe.subset <- spe[, colData(spe)$sample_id == id]
    p <- make_escheR(spe.subset) |>
        add_fill(var="clust_HARMONY_M0_lam0.4_k50_res0.4", point_size=1.5)+
        scale_fill_manual(values = pal)
    print(p)
}
dev.off()


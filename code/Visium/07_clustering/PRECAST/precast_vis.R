setwd('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("PRECAST")
    library("dplyr")
    library("purrr")
    library("tidyverse")
    library("spatialLIBD")
    library("gridExtra")
    library("ggspavis")
    library("escheR")
})

load(here("processed-data","Visium", "06_batch_correction", "spe_harmony.Rdata"))
colnames(spe) <- spe$key

# drop duplicated colData, if any
duplicated_cols <- duplicated(colnames(colData(spe)))
filtered_colData <- colData(spe)[, !duplicated_cols]
colData(spe) <- filtered_colData

K <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
load(file = here("processed-data", "Visium", "07_clustering", "PRECAST","SVGs", paste0("PRECASTObj_",K,".Rdata")))

PRECASTObj <- SelectModel(PRECASTObj)
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

# Merge with spe object
cluster_df <- seuInt@meta.data |>
    mutate(cluster = factor(cluster)) |>
    rename_with(~ paste0("PRECAST_", .x)) |>
    rownames_to_column(var = "key")

col_data_df <- colData(spe) |>
    data.frame() |>
    left_join(cluster_df, by="key")

rownames(col_data_df) <- colnames(spe)
colData(spe)$PRECAST_cluster <- col_data_df$PRECAST_cluster

precast_name <- paste0("PRECAST_clusters_", K)

cluster_export(
    spe,
    "PRECAST_cluster",
    cluster_dir = here::here("processed-data", "Visium", "07_clustering", "PRECAST","SVGs", "cluster_csv", precast_name)
)

# loop through each sample_id and create a pdf on a separete page
pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(spe.subset$PRECAST_cluster)))
pdf(file = here::here("plots", "Visium", "07_clustering", "PRECAST", "SVGs", paste0(precast_name, ".pdf")), width = 10, height = 10)
for (id in unique(colData(spe)$sample_id)) {
    spe.subset <- spe[, colData(spe)$sample_id == id]
    p <- make_escheR(spe.subset) |>
        add_fill(var="PRECAST_cluster")+
        scale_fill_manual(values = pal)
    print(p)
}
dev.off()
setwd('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("PRECAST")
    library("tictoc")
    library("dplyr")
    library("purrr")
    library("tidyverse")
    library("spatialLIBD")
    library("gridExtra")
    library("ggspavis")
})

load(here("processed-data", "03_qc_metrics", "spe_discarded.Rdata"))
colnames(spe) <- spe$key

# drop duplicated colData, if any
duplicated_cols <- duplicated(colnames(colData(spe)))
filtered_colData <- colData(spe)[, !duplicated_cols]
colData(spe) <- filtered_colData

K <- as.numeric(Sys.getenv("SGE_TASK_ID"))
load(file = here("processed-data", "04_dimred_clustering", "PRECAST", paste0("PRECASTObj_",K,".Rdata")))

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
    cluster_dir = here::here("processed-data", "04_dimred_clustering", "PRECAST", precast_name)
)


pdf(file = here::here("plots", "04_dimred_clustering", "PRECAST", paste0(precast_name, ".pdf")), width = 21, height = 12)

p1 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[1], clustervar = "PRECAST_cluster",point_size = 1.5)
p2 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[2], clustervar = "PRECAST_cluster",point_size = 1.5)
p3 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[3], clustervar = "PRECAST_cluster",point_size = 1.5)
p4 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[4], clustervar = "PRECAST_cluster",point_size = 1.5)
p5 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[5], clustervar = "PRECAST_cluster",point_size = 1.5)
p6 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[6], clustervar = "PRECAST_cluster",point_size = 1.5)
p7 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[7], clustervar = "PRECAST_cluster",point_size = 1.5)
p8 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[8], clustervar = "PRECAST_cluster",point_size = 1.5)

(p1|p2|p3|p4)/(p5|p6|p7|p8)
dev.off()
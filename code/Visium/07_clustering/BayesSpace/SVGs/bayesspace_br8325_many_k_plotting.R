suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("RColorBrewer")
    library("ggplot2")
    library("patchwork")
})

load(here("processed-data","Visium", "06_batch_correction", "spe_harmony.Rdata"))
dim(spe)

#subset to brnum 9280
spe <- spe[, colData(spe)$sample_id == "Br8325"]
spe

# get folders in cluster_Csv
bs_folders <- list.files(here::here("processed-data","Visium", "07_clustering", "BayesSpace","SVGs","Br8325"), full.names = TRUE)
bs_folders
# [1] "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/08_clustering/BayesSpace/SVGs/cluster_csv/BayesSpace_10"
# [2] "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/08_clustering/BayesSpace/SVGs/cluster_csv/BayesSpace_12"

# get the number of clusters at the end of the folder name
bs_k <- gsub(".*BayesSpace_", "", bs_folders)

# loop through each folder, open the csv inside, and add clusters to spe colData
for (i in seq_along(bs_folders)) {
    bs_folder <- bs_folders[i]
    bs_csv <- list.files(bs_folder, full.names = TRUE)
    bs_csv <- bs_csv[grepl("csv", bs_csv)]
    bs_csv <- bs_csv[1]
    bs_df <- read.csv(bs_csv)

    colData(spe)[[paste0("BS_k", bs_k[i])]] <- factor(bs_df$cluster)
}

colnames(colData(spe))

# loop through each "BS_i" column and create a pdf on a separete page for sample_id. One pdf per BS_i
library("escheR")
for (i in seq_along(bs_folders)) {
    clustV <- paste0("BS_k", bs_k[i])
    pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(colData(spe)[[clustV]])))
    pdf(file = here::here("plots", "Visium", "07_clustering", "BayesSpace","SVGs", "Br8325", paste0("BS_k", bs_k[i], "_stitched.pdf")), width = 10, height = 10)
    for (sample_id in unique(colData(spe)$sample_id)) {
        spe.subset <- spe[, colData(spe)$sample_id == sample_id]
        p <- make_escheR(spe.subset) |>
            add_fill(var=clustV, point_size=1.75) +
            scale_fill_manual(values = pal)
        print(p)
    }
    dev.off()
}


# drop spe$exclude_overlapping and NAs
spe.test <- spe[, !spe$exclude_overlapping & !is.na(spe$exclude_overlapping)]

for (i in seq_along(bs_folders)) {
    clustV <- paste0("BS_k", bs_k[i])
    pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(colData(spe.test)[[clustV]])))
    pdf(file = here::here("plots", "Visium", "07_clustering", "BayesSpace", "Br8325", paste0("LIBD_BS_k", bs_k[i], "_stitched.pdf")), width = 10, height = 10)
    p <- vis_clus(
        spe.test,
        clustervar=clustV,
        colors = pal,
        point_size = 3,
        is_stitched = TRUE,
        )
        print(p)
    dev.off()
}
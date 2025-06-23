setwd('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/')
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
})

load(here("processed-data","Visium", "06_batch_correction", "spe_harmony.Rdata"))
dim(spe)
 

#subset to brnum 9280
spe <- spe[, colData(spe)$sample_id == "Br8325"]

k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

spe <- spatialPreprocess(spe, platform="Visium", skip.PCA=TRUE)

colData(spe)$row <- spe$array_row
colData(spe)$col <- spe$array_col

metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

message("Running spatialCluster()")
Sys.time()
set.seed(2)
spe <- spatialCluster(spe, use.dimred = "PCA-HARMONY_sample", q = k, nrep=10000, burn.in=100, gamma=3)
Sys.time()

bayesSpace_name <- paste0("BayesSpace_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = here::here("processed-data","Visium", "07_clustering", "BayesSpace","SVGs","Br8325")
)

clustV <- bayesSpace_name


# loop through each sample_id and create a pdf on a separete page
pdf(file = here::here("plots","Visium", "07_clustering", "BayesSpace", "SVGs", "Br8325", paste0(bayesSpace_name, ".pdf")), width = 10, height = 10)
for (sample_id in unique(colData(spe)$sample_id)) {
    spe.subset <- spe[, colData(spe)$sample_id == sample_id]
    p <- make_escheR(spe.subset) |>
        add_fill(var=clustV)
    print(p)
}
dev.off()

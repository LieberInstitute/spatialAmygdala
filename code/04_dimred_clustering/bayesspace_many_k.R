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
    library("Polychrome")
})

load(here("processed-data", "03_qc_metrics", "spe_discarded.Rdata"))
dim(spe)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))


spe <- spatialPreprocess(spe, platform="Visium", n.PCs=7, n.HVGs=2000, log.normalize=TRUE)

colData(spe)$row <- spe$array_row
colData(spe)$col <- spe$array_col

metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

message("Running spatialCluster()")
Sys.time()
set.seed(2)
spe <- spatialCluster(spe, q = k, nrep=10000, burn.in=100)
Sys.time()

bayesSpace_name <- paste0("BayesSpace_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = here::here("processed-data", "04_dimred_clustering", "BayesSpace", bayesSpace_name)
)

cols <- Polychrome::palette36.colors(k)
names(cols) <- sort(unique(spe$spatial.cluster))
clustV <- bayesSpace_name

pdf(file = here::here("plots", "04_dimred_clustering", "BayesSpace", paste0(bayesSpace_name, ".pdf")), width = 21, height = 20)

p1 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[1], clustervar = clustV, colors = cols, point_size = 5)
p2 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[2], clustervar = clustV, colors = cols, point_size = 5)
p3 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[3], clustervar = clustV, colors = cols, point_size = 5)
p4 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[4], clustervar = clustV, colors = cols, point_size = 5)
p5 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[5], clustervar = clustV, colors = cols, point_size = 5)
p6 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[6], clustervar = clustV, colors = cols, point_size = 5)
p7 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[7], clustervar = clustV, colors = cols, point_size = 5)
p8 <- vis_clus(spe = spe, sampleid = unique(colData(spe)$sample_id)[8], clustervar = clustV, colors = cols, point_size = 5)

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2)
dev.off()

library("SpatialExperiment")
library("here")
library("scater")
library("nnSVG")
library("scran")
library("BayesSpace")

# save directiories
plot_dir = here("plots", "04_dimred_clustering")
processed_dir = here("processed-data", "04_dimred_clustering")

# load object
load(here("processed-data","03_qc_metrics","spe_discarded.Rdata"))
spe

# --------------- preprocess ------------------

set.seed(102)
spe <- spatialPreprocess(spe, platform="Visium")
spe

# add 'row' and 'col' colData for spatialClustering
colData(spe)$row <- colData(spe)$array_row
colData(spe)$col <- colData(spe)$array_col

spe <- spatialCluster(spe, q=10, nrep=100, gamma=3, save.chain=TRUE, burn.in = 10)
spe

# plot
pdf(here(plot_dir,"BS_test_clusters_k10.pdf"))
clusterPlot(spe)
dev.off()

# save
save(spe, file=here(processed_dir,"spe_clusters_k10.Rdata"))
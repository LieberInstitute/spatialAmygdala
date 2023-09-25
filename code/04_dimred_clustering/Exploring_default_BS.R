library("SpatialExperiment")
library("here")
library("scater")
library("nnSVG")
library("scran")
library("BayesSpace")
library("spatialLIBD")
library("patchwork")

# save directiories
plot_dir = here("plots", "04_dimred_clustering")
processed_dir = here("processed-data", "04_dimred_clustering")

# load object
load(here("processed-data","04_dimred_clustering","spe_clusters_k10.Rdata"))
spe

# Spotplots
pdf(width=20, height=10, here(plot_dir,"BS_387-C1_spotplot.pdf"))
vis_clus(
    spe = spe,
    clustervar = "spatial.cluster",
    sampleid = "V13M06-387_C1"
  )
dev.off()

pdf(width=20, height=10, here(plot_dir,"BS_387-D1_spotplot.pdf"))
vis_clus(
    spe = spe,
    clustervar = "spatial.cluster",
    sampleid = "V13M06-387_D1"
  )
dev.off()
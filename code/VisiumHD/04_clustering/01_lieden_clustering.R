library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialFeatureExperiment)
library(scater)
library(scran)
library(Voyager)
library(ggspavis)
library(igraph)
library(scry)
library(harmony)

# set directories
plots_dir <- here("plots", "VisiumHD", "04_clustering")
processed_dir <- here("processed-data", "VisiumHD", "03_quality_control")

# load
sfe <- readRDS(here(processed_dir, "sfe_qc.rds"))

# ========== Lieden clustering ==========

set.seed(1234)

# add donor based on Sample
sfe$donor <- as.factor(ifelse(grepl("Br8325", sfe$Sample), "Br8325", "Br9280"))

# normalization
sfe <- logNormCounts(sfe)

# HVG
dec <- modelGeneVar(sfe)
top.hvgs <- getTopHVGs(dec, n=2000)

# PCA
sfe <- runPCA(sfe, subset_row=top.hvgs, ncomponents=50)

# Harmony
sfe <- RunHarmony(sfe, group.by.vars="donor", dimred="PCA")

# clustering
g <- buildSNNGraph(sfe, use.dimred="HARMONY", type="jaccard")
k <- cluster_leiden(g, objective_function="modularity", resolution=0.7)

# table of clusters
table(k$membership)
#     1     2     3     4     5     6     7     8     9    10 
# 24369   911 84875 36469 58486  5061  5984 24010 37589  7232 

sfe$leiden <- as.factor(k$membership)

# UMAP
sfe <- runUMAP(sfe, dimred="HARMONY")

#  ============= Plotting =============


pdf(here(plots_dir, "lieden_clustering.pdf"), width=10, height=5)
p1 <- plotUMAP(sfe, colour_by="leiden") +
  theme(legend.position = "right") +
  guides(colour=guide_legend(override.aes=list(size=5))) +
  ggtitle("Leiden Clustering (res=0.7)") +
  theme(plot.title = element_text(hjust = 0.5, size=14))
p2 <- plotUMAP(sfe, colour_by="Sample") +
    theme(legend.position = "right") +
    guides(colour=guide_legend(override.aes=list(size=5))) +
    ggtitle("Leiden Clustering (res=0.7)") +
    theme(plot.title = element_text(hjust = 0.5, size=14))
p1 + p2
dev.off()

# spot plots
pdf(here(plots_dir, "lieden_clustering_spotplots.pdf"), width=9, height=6)
plotSpatialFeature(sfe, features = "leiden", colGeometryName = "cellSeg") + 
    ggtitle("Leiden Clustering") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.key.width = grid::unit(0.5, "lines"), 
          legend.key.height = grid::unit(1, "lines"))
dev.off()

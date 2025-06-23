library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("SingleR")

# save directories
processed_dir <- here("processed-data", "Xenium", "04_clusterig")
plot_dir <- here("plots", "Xenium", "04_clustering")


# load xenium data
load(here("processed-data","Xenium", "03_quality_control", "spe_normcounts.Rdata"))
spe

# load the snRNA-seq data
sce <- readRDS(here("processed-data", "snRNAseq", "sce.human_all_genes.rds"))
sce

# ======= Label transfer =======

# set rownames to symbols
rownames(spe) <- rowData(spe)$Symbol

# subset to genes in both datasets
common_genes <- intersect(rownames(spe), rownames(sce))
spe <- spe[common_genes, ]
sce <- sce[common_genes, ]

# re-normalize the snRNAseq data
spe <- computeLibraryFactors(spe)
spe <- spe[, sizeFactors(spe) > 0]
spe <- logNormCounts(spe)

sce <- logNormCounts(sce)

dim(spe)
#[1]     366 1017874

dim(sce)
#[1]   366 15511

# Broad cell type label transfer ===
pred.broad <- SingleR(test=spe, ref=sce, labels=sce$broad_celltype, de.method="wilcox")
table(pred.broad$labels)

# Fine cell type label transfer ===
pred.fine <- SingleR(test=spe, ref=sce, labels=sce$fine_celltype, de.method="wilcox")
table(pred.fine$labels)

# save to csv
write.csv(pred.broad, file=here("processed-data", "Xenium", "06_label_transfer", "pred_broad_celltype.csv"), row.names=FALSE)
write.csv(pred.fine, file=here("processed-data", "Xenium", "06_label_transfer", "pred_fine_celltype.csv"), row.names=FALSE)
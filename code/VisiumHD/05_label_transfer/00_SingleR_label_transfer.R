library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("SingleR")
library("BiocParallel")

plot_dir <- here("plots","VisiumHD","05_label_transfer", "segmentations")


spe <- readRDS(here("processed-data/VisiumHD/02_build_spe/Br9280_CeA_Spatial.Polygons_spe.rds"))
spe
# dim: 18085 50908 
# metadata(0):
# assays(1): counts
# rownames(18085): SAMD11 NOC2L ... MT-ND6 MT-CYB
# rowData names(1): gene_id
# colnames(50908): cellid_000000004-1 cellid_000000007-1 ...
#   cellid_000060097-1 cellid_000060106-1
# colData names(9): orig.ident nCount_Spatial.008um ... sample_id
#   in_tissue
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x y
# imgData names(4): sample_id image_id data scaleFactor


# load the snRNA-seq data
sce <- readRDS(here("processed-data", "snRNAseq", "sce.human_all_genes.rds"))
sce
# class: SingleCellExperiment 
# dim: 36601 15511 
# metadata(0):
# assays(1): counts
# rownames(36601): A1BG A1BG-AS1 ... ZYX ZZEF1
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(15511): cell_2 cell_4 ... cell_16967 cell_16969
# colData names(33): batch orig.ident ... fine_celltype ident
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):


# ======= Label transfer ========

# normalize counts 
spe <- logNormCounts(spe)
sce <- logNormCounts(sce)

# subset to genes in both datasets
common_genes <- intersect(rownames(spe), rownames(sce))
spe <- spe[common_genes, ]
sce <- sce[common_genes, ]

dim(spe)
# [1] 18076 50908

dim(sce)
# [1] 18076 15511

# Broad cell type label transfer ===
pred.broad <- SingleR(test=logcounts(spe), ref=logcounts(sce), labels=sce$broad_celltype, de.method="wilcox",
    BPPARAM=MulticoreParam(20), de.n=100)
table(pred.broad$labels)

# Fine cell type label transfer ===
pred.fine <- SingleR(test=logcounts(spe), ref=logcounts(sce), labels=sce$fine_celltype, de.method="wilcox",
    BPPARAM=MulticoreParam(20), de.n=100)
table(pred.fine$labels)

# save to csv
write.csv(pred.broad, file=here("processed-data", "VisiumHD", "05_label_transfer", "SingleR", "HDsegmentations_pred_broad_celltype.csv"), row.names=FALSE)
write.csv(pred.fine, file=here("processed-data", "VisiumHD", "05_label_transfer", "SingleR", "HDsegmentations_pred_fine_celltype.csv"), row.names=FALSE)
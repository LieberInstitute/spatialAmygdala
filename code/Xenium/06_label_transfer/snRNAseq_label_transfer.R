library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("SingleR")
library("BiocParallel")


# load xenium data
spe <- readRDS(here("processed-data","Xenium", "04_dim_reduction", "spe_xenium_5um_harmonized_singlecell.rds"))
spe

# load the snRNA-seq data
sce <- readRDS(here("processed-data", "snRNAseq", "sce.human_all_genes.rds"))
sce

logcounts(sce) <- logNormCounts(sce)
logcounts(spe) <- assay(spe, "cell_normcounts")

class(logcounts(sce))
class(logcounts(spe))

# ======= Label transfer =======

# subset to genes in both datasets
common_genes <- intersect(rownames(spe), rownames(sce))
spe <- spe[common_genes, ]
sce <- sce[common_genes, ]

dim(spe)
#[1]     366 1017874

dim(sce)
#[1]   366 15511

# Broad cell type label transfer ===
pred.broad <- SingleR(test=logcounts(spe), ref=logcounts(sce), labels=sce$broad_celltype, de.method="wilcox",
    BPPARAM=MulticoreParam(20), de.n=100)
table(pred.broad$labels)

# Fine cell type label transfer ===
pred.fine <- SingleR(test=logcounts(spe), ref=logcounts(sce), labels=sce$fine_celltype, de.method="wilcox",
    BPPARAM=MulticoreParam(20), de.n=100)
table(pred.fine$labels)

# save to csv
write.csv(pred.broad, file=here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_broad_celltype_n100.csv"), row.names=FALSE)
write.csv(pred.fine, file=here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_fine_celltype_n100.csv"), row.names=FALSE)
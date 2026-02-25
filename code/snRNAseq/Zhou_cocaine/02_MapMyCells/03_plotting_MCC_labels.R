library(here)
library(SingleCellExperiment)
library(scater)
library(scran)

# load zhou RDS
sce <- readRDS(here("processed-data","snRNAseq","Zhou_cocaine","Zhou_sce_NMF50_with_MapMyCells.rds"))
sce
# class: SingleCellExperiment 
# dim: 17297 163003 
# metadata(1): NMF_basis
# assays(2): counts logcounts
# rownames(17297): AABR07000156.1 Lrp11 ... AABR07043200.1 Pomp
# rowData names(0):
# colnames(163003): AAACCCAAGAAACCCG-1_1 AAACCCACAAAGCACG-1_1 ...
#   TTTGTTGTCTTCGTAT-1_19 TTTGTTGTCTTCTGGC-1_19
# colData names(49): orig.ident nCount_RNA ... mmc_cluster_alias
#   mmc_cluster_bootstrapping_probability
# reducedDimNames(1): NMF
# mainExpName: RNA
# altExpNames(2): SCT integrated


# ======= Generste PC and UMAP spaces =======

set.seed(1234)
# get 3000 HVGs
dec <- modelGeneVar(sce)
top_hvgs <- getTopHVGs(dec, n=3000)

# run PCA
sce <- runPCA(sce, subset_row=top_hvgs, ncomponents=75)

# run Harmony
library(harmony)
sce <- RunHarmony(sce, "sample", assay.use="logcounts", reduction="PCA", reduction.save="HARMONY")

# run UMAP on Harmony embeddings
sce <- runUMAP(sce, dimred="HARMONY", min_dist=0.3)


# ======= UMAP plots with MapMyCells labels =======

pdf(here("plots","snRNAseq","Zhou_cocaine","Zhou_sce_NMF50_MapMyCells_UMAP_by_mmc_cluster_alias.pdf"), width=8, height=6)
plotUMAP(sce, colour_by="mmc_class_name") +
  ggtitle("MapMyCells cluster alias")
dev.off()

table(sce$mmc_subclass_name)

# get number of subclasses less than 500 cells
subclass_counts <- table(sce$mmc_subclass_name)
small_subclasses <- names(subclass_counts[subclass_counts < 500])
length(small_subclasses)

large_subclasses <- names(subclass_counts[subclass_counts >= 500])
length(large_subclasses)

sce.subset <- sce 

# drop small subclasses
sce.subset <- sce.subset[, !(sce.subset$mmc_subclass_name %in% small_subclasses)]


pdf(here("plots","snRNAseq","Zhou_cocaine","Zhou_sce_NMF50_MapMyCells_UMAP_by_mmc_subclass_name_filtered.pdf"), width=10, height=8)
plotUMAP(sce.subset, colour_by="mmc_subclass_name") +
    ggtitle("MapMyCells subclass name (filtered)")
dev.off()

# plot with a color palette for 40 different cell types using polychrome (A LOT)
library(RColorBrewer)
n <- length(large_subclasses)
getPalette <- colorRampPalette(brewer.pal(8, "Paired"))
colors <- getPalette(n) 

pdf(here("plots","snRNAseq","Zhou_cocaine","Zhou_sce_NMF50_MapMyCells_UMAP_by_mmc_subclass_name_filtered_colorblind.pdf"), width=10, height=8)
plotUMAP(sce.subset, colour_by="mmc_subclass_name") +
    scale_color_manual(values = colors) +
    ggtitle("MapMyCells subclass name (filtered, colorblind)")
dev.off()


# plot NMF components 15, 17, 3
library(NMFscape)
p1 <- vizUMAP(sce, program = 3) +
    ggtitle("NMF Component 3")
p2 <- vizUMAP(sce, program = 15) +
    ggtitle("NMF Component 15")
p3 <- vizUMAP(sce, program = 17) +
    ggtitle("NMF Component 17")

pdf(here("plots","snRNAseq","Zhou_cocaine","Zhou_sce_NMF50_UMAP_NMF_components_3_15_17.pdf"), width=12, height=4)
library(cowplot)
plot_grid(p1, p2, p3, nrow=1)
dev.off()


# dot plot across MapMyCells subclasses
pdf(here("plots","snRNAseq","Zhou_cocaine","Zhou_sce_NMF50_MapMyCells_dotplot_NMF_components_filtered.pdf"), width=12, height=20)
plotProgramDots(sce.subset, group = "mmc_subclass_name") 
dev.off()


# DE analysis
deps <- FindAllDEPs(sce.subset, cell_type_col = "mmc_subclass_name")

# Heatmap of log2 fold change with stars for significant enrichments
pdf(here("plots","snRNAseq","Zhou_cocaine","Zhou_sce_NMF50_MapMyCells_DEPs_heatmap_filtered.pdf"), width=10, height=18)
plotDEPsHeatmap(deps, logfc_threshold = 2, star_size=14)
dev.off()
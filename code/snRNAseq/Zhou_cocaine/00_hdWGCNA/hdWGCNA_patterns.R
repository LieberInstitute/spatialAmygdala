library("here")
library("Seurat")
library("SingleCellExperiment")
library("hdWGCNA")


# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 20)

# load data
rat.amy <- readRDS(here("processed-data","snRNAseq", "GSE212415_seurat.rds"))
rat.amy

# add cell type information
rat.amy$idents <- Idents(rat.amy)

# = setup for hdWGCNA = #
seurat_obj <- SetupForWGCNA(
    rat.amy,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "Cocaine" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("idents", "label", "sample"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'pca', # select the dimensionality reduction to perform KNN on
    k = 25, # nearest-neighbors parameter
    max_shared = 10, # maximum number of shared cells between two metacells
    ident.group = 'idents' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
    seurat_obj,
    layer = "data", # use the normalized data slot
    assay = "RNA" # specify the assay to pull data from
)


# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
pdf( file = here("plots","snRNAseq","Zhou_cocaine","hdWGCNA_softpower_plots.pdf"), width = 8, height = 6)
wrap_plots(plot_list, ncol=2)
dev.off()
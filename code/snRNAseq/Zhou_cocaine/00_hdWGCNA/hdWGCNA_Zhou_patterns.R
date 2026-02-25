library("here")
library("Seurat")
library("SingleCellExperiment")
library("hdWGCNA")
library("patchwork")
library("dplyr")
library("UCell")
library("ggplot2")


# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 30)

# load data
rat.amy <- readRDS(here("processed-data","snRNAseq", "GSE212415_seurat.rds"))
rat.amy

# add cell type information
rat.amy$idents <- Idents(rat.amy)
#  [1] InhNeuron        Astrocytes       ExNeuron         Nos1+           
#  [5] Microglia        Oligodendrocytes Chat+            Cck+/Vip+       
#  [9] OPC              Sst+             Reln+            Endothelial     
# [13] Pvalb+          
# 13 Levels: Cck+/Vip+ Astrocytes Oligodendrocytes InhNeuron OPC ... Pvalb+

# new column broad_celltype based on ExNeuron, InhNeuron, and Glia. Where InhNeuron includes all inhibitory neuron subtypes, and Glia includes all glial cell types.
rat.amy$broad_celltype <- case_when(
  rat.amy$idents == "ExNeuron" ~ "Excitatory Neuron",
  rat.amy$idents %in% c("InhNeuron", "Nos1+", "Cck+/Vip+", "Sst+", "Reln+", "Pvalb+", "Chat+") ~ "Inhibitory Neuron",
  rat.amy$idents %in% c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC", "Endothelial") ~ "Glia",
  TRUE ~ "Other"
)


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
    group.by = c("broad_celltype", "label", "sample"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'pca', # select the dimensionality reduction to perform KNN on
    k = 50, # nearest-neighbors parameter
    max_shared = 10, # maximum number of shared cells between two metacells
    ident.group = 'broad_celltype' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = "Inhibitory Neuron", # the name of the group of interest in the group.by column
    group.by='broad_celltype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
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
pdf( file = here("plots","snRNAseq","Zhou_cocaine","INH_hdWGCNA_softpower_plots.pdf"), width = 8, height = 6)
wrap_plots(plot_list, ncol=2)
dev.off()


# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'Zhou', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE,
  deepSplit = 3,         # higher = more modules (try 2–4)
  minModuleSize = 15,    # lower = more modules (try 15–30)
  mergeCutHeight = 0.05  # lower = less merging (try 0.05–0.20)
)

pdf( file = here("plots","snRNAseq","Zhou_cocaine","INH_hdWGCNA_dendrogram.pdf"), width = 8, height = 6)
PlotDendrogram(seurat_obj, main='Zhou hdWGCNA Dendrogram')
dev.off()


# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="sample"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(seurat_obj)

# rename the modules
seurat_obj <- ResetModuleNames(seurat_obj, new_name = "Zhou-M")

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

pdf( file = here("plots","snRNAseq","Zhou_cocaine","INH_hdWGCNA_kME_plots.pdf"), width = 12, height = 8)
p
dev.off()


# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])

saveRDS(seurat_obj, file=here("processed-data","snRNAseq","Zhou_cocaine","INH_hdWGCNA_Zhou_object.rds"))


# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)


# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
pdf( file = here("plots","snRNAseq","Zhou_cocaine","INH_hdWGCNA_hME_featureplots.pdf"), width = 16, height = 12)
wrap_plots(plot_list, ncol=6)
dev.off()

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
pdf( file = here("plots","snRNAseq","Zhou_cocaine","INH_hdWGCNA_hubscore_featureplots.pdf"), width = 16, height = 12)
wrap_plots(plot_list, ncol=6)
dev.off()


# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)


# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'idents')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
pdf( file = here("plots","snRNAseq","Zhou_cocaine","INH_hdWGCNA_hME_dotplot.pdf"), width = 8, height = 6)
p
dev.off()
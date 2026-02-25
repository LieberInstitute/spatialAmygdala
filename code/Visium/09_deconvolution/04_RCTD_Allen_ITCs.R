library(SingleCellExperiment)
library(here)
library(spacexr)

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")

# load
sce.libd <- readRDS(here("processed-data", "snRNAseq", "sce_amy_macaque.rds"))
sce.libd

rownames(sce.libd) <- rowData(sce.libd)$gene_id

# drop TSHZ1 cell from fine_celltype using grep
sce.libd <- sce.libd[, !grepl("TSHZ1", sce.libd$fine_celltype)]


# --- Load SCEs ---
sce.itc <- readRDS(here("processed-data/snRNAseq/BICCN/WHB_amygdala_Final_ITC_subsets.rds"))
sce.itc
# class: SingleCellExperiment 
# dim: 59357 66072 
# metadata(0):
# assays(1): counts
# rownames(59357): ENSG00000000003 ENSG00000000005 ... ENSG00000288643
#   ENSG00000288645
# rowData names(1): gene_symbol
# colnames(66072): 10X386_2:CCCTCAAAGTCCCAAT 10X383_5:ATGCCTCCATCCTATT
#   ... 10X218_5:CGACAGCAGATTAGTG 10X204_3:TTTCACATCGACGACC
# colData names(30): cell_barcode library_label ...
#   neurotransmitter_color region_of_interest_color
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):

# Access counts matrix (assuming it's called "counts" assay)
counts_mat <- assay(sce.itc, "counts")

# Check if all entries are integers
all_integers <- all(counts_mat == round(counts_mat))

table(sce.itc$supercluster)
    #       Amygdala excitatory         Cerebellar inhibitory 
    #                        48                           121 
    #           CGE interneuron Eccentric medium spiny neuron 
    #                      1894                         28357 
    # LAMP5-LHX6 and Chandelier             Lower rhombic lip 
    #                       383                           544 
    #           Mammillary body           Medium spiny neuron 
    #                       155                            52 
    #           MGE interneuron   Midbrain-derived inhibitory 
    #                      1734                          1050 
    #             Miscellaneous                      Splatter 
    #                         1                          2713 
    #       Thalamic excitatory             Upper rhombic lip 
    #                       731                          1493 


# keep only Eccentric medium spiny neuron
sce.itc <- sce.itc[, sce.itc$supercluster == "Eccentric medium spiny neuron"]
sce.itc
# class: SingleCellExperiment 
# dim: 59357 28357 
# metadata(0):
# assays(2): counts logcounts
# rownames(59357): TSPAN6 TNMD ... AC114982.3 AC084756.2
# rowData names(2): gene_symbol gene_id
# colnames(28357): 10X356_8:ATTACCTAGGAGCAAA 10X356_8:ATCGCCTGTTGGACTT
#   ... 10X349_4:GACTCAACAGACCGCT 10X354_7:CGCATGGGTCAGTCGC
# colData names(32): cell_barcode library_label ... sizeFactor
#   louvain_cluster
# reducedDimNames(4): PCA UMAP harmony UMAP_harmony
# mainExpName: NULL
# altExpNames(0):


# change libd to gene symbols
rownames(sce.libd) <- rowData(sce.libd)$gene_name

# subset sce.itc and sce.libd to common genes
common_genes <- intersect(rownames(sce.itc), rownames(sce.libd))
sce.itc <- sce.itc[common_genes, ]
sce.libd <- sce.libd[common_genes, ]

# add $fine_celltype from to sce.itc
sce.itc$fine_celltype <- sce.itc$cluster



# keep only common columns in colData prior to merge
common_cols <- intersect(colnames(colData(sce.libd)), colnames(colData(sce.itc)))
common_cols

# subset colData to common colData. make sure its colData
colData(sce.libd) <- colData(sce.libd)[, common_cols]
colData(sce.itc) <- colData(sce.itc)[, common_cols]

# set gene_id to NULL from rowData prior to merge
rowData(sce.libd)$gene_id <- NULL
rowData(sce.itc)$gene_id <- NULL

# drop reducedDims
reducedDims(sce.libd) <- NULL
reducedDims(sce.itc) <- NULL

# merge sce objects
sce.combined <- cbind(sce.libd, sce.itc)

unique(sce.combined$fine_celltype)
#  [1] "EMSN_426" "EMSN_227" "EMSN_226" "EMSN_225" "EMSN_230" "EMSN_231"
#  [7] "EMSN_232" "EMSN_233" "EMSN_234" "EMSN_228" "EMSN_229" "EMSN_222"
# [13] "EMSN_223" "EMSN_224"

# load
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
spe <- spe[!duplicated(rownames(spe)), ]

rownames(spe) <- rowData(spe)$gene_name

# ========== RCTD ==========
rctd_data <- createRctd(spe, sce.combined, cell_type_col="fine_celltype")
res <- runRctd(rctd_data, max_cores=30, rctd_mode="multi", max_multi_types=5)

# save
saveRDS(res, here(processed_dir, "rctd_Allen_ITCs_results.rds"))
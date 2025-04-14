suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SingleCellExperiment")
    library("spatialLIBD")
    library("RColorBrewer")
    library("ggplot2")
    library("patchwork")
    library("MetaNeighbor")
})

# === set up directories ===
processed_dir <- here("processed-data", "snRNAseq")
plot_dir <- here("plots", "Visium", "09_marker_genes", "K99")

sce.totty <- readRDS(here(processed_dir, "sce.human_all_genes.rds"))
sce.totty
# class: SingleCellExperiment 
# dim: 21844 20625 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(21844): AL627309.1 AL627309.5 ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(20625): AAACCCAAGCTAAATG-1 AAACCCACAGGTCCCA-1 ...
#   TTTGTTGTCGGACTGC-1 TTTGTTGTCGTTGTTT-1
# colData names(34): orig.ident nCount_originalexp ... broad_celltype
#   clusters
# reducedDimNames(2): PCA UMAP
# mainExpName: NULL
# altExpNames(0):

sce.zhuo <- readRDS(here(processed_dir, "Zhuo_cocaine_rats_excitatory_clustered.rds"))
sce.zhuo
# class: SingleCellExperiment 
# dim: 17297 23943 
# metadata(0):
# assays(2): counts logcounts
# rownames(17297): AABR07000156.1 Lrp11 ... AABR07043200.1 Pomp
# rowData names(0):
# colnames(23943): AAACCCAGTGCATACT-1_1 AAACCCAGTGTCTAAC-1_1 ...
#   TTTGTTGGTTGCTCAA-1_19 TTTGTTGTCTTCTGGC-1_19
# colData names(40): orig.ident nCount_RNA ... lieden_k25 lieden_k30
# reducedDimNames(3): PCA UMAP GLM-PCA
# mainExpName: RNA
# altExpNames(0):



# =========================================
# subset to 1-1 orthologs
# =========================================


 # ====== Subset to only 1:1 orthologs ======

 method <- "gprofiler"

rat_orthos <- orthogene::convert_orthologs(gene_df = counts(sce.zhuo),
                                         gene_input = "rownames", 
                                         gene_output = "rownames", 
                                         input_species = "rat",
                                         output_species = "human",
                                         non121_strategy = "drop_both_species",
                                         method = method) 

# =========== REPORT SUMMARY ===========
# Total genes dropped after convert_orthologs :
#    4,041 / 17,297 (23%)
# Total genes remaining after convert_orthologs :
#    13,256 / 17,297 (77%)


# Extract rat gene names that have 1:1 orthologs
ortholog_genes <- rownames(rat_orthos)

# Subset sce.zhuo to only these genes
sce.zhuo.orthos <- sce.zhuo[names(ortholog_genes), ]

rowData(sce.zhuo.orthos)$ortholog <- ortholog_genes

rownames(sce.zhuo.orthos) <- rowData(sce.zhuo.orthos)$ortholog


 # get only orthologs that are in the human dataset
 valid.human.orthologs <- rownames(rat_orthos)[rownames(rat_orthos) %in% rownames(sce.totty)]
 length(valid.human.orthologs)
# [1] 12069

# subset to only 1:1 orthologs
sce.human.ortho <- sce.totty[valid.human.orthologs ,]
sce.human.ortho
# class: SingleCellExperiment 
# dim: 12069 23943 
# metadata(0):
# assays(2): counts logcounts
# rownames(12069): LRP11 PCMT1 ... COL23A1 POMP
# rowData names(1): ortholog
# colnames(23943): AAACCCAGTGCATACT-1_1 AAACCCAGTGTCTAAC-1_1 ...
#   TTTGTTGGTTGCTCAA-1_19 TTTGTTGTCTTCTGGC-1_19
# colData names(40): orig.ident nCount_RNA ... lieden_k25 lieden_k30
# reducedDimNames(3): PCA UMAP GLM-PCA
# mainExpName: RNA
# altExpNames(0):

sce.zhuo.ortho <- sce.zhuo.orthos[valid.human.orthologs,]
sce.zhuo.ortho
# class: SingleCellExperiment 
# dim: 12069 23943 
# metadata(0):
# assays(2): counts logcounts
# rownames(12069): LRP11 PCMT1 ... COL23A1 POMP
# rowData names(1): ortholog
# colnames(23943): AAACCCAGTGCATACT-1_1 AAACCCAGTGTCTAAC-1_1 ...
#   TTTGTTGGTTGCTCAA-1_19 TTTGTTGTCTTCTGGC-1_19
# colData names(40): orig.ident nCount_RNA ... lieden_k25 lieden_k30
# reducedDimNames(3): PCA UMAP GLM-PCA
# mainExpName: RNA
# altExpNames(0):


sce.totty.ortho <- sce.human.ortho

# subset to human "Excitatory" neurons in $broad_celltype
sce.totty.ortho <- sce.totty.ortho[, sce.totty.ortho$broad_celltype == "Excitatory"]

# refactor fine_celltype
sce.totty.ortho$fine_celltype <- factor(sce.totty.ortho$fine_celltype)

#==========================================
# Combined the datasets
#==========================================

# common naming for cell type clusters
sce.totty.ortho$clusters <- sce.totty.ortho$fine_celltype
sce.zhuo.ortho$clusters <- sce.zhuo.ortho$lieden_k25

# Assign rowRanges from sce.totty to sce.zhuo
rowRanges(sce.zhuo.ortho) <- rowRanges(sce.totty.ortho)

# Add "Rat_" and "Human_" prefixes to the clusters
#sce.totty.ortho$clusters <- paste0("Human_", sce.totty.ortho$clusters)
#sce.zhuo.ortho$clusters <- paste0("Rat_", sce.zhuo.ortho$clusters)

# add a species column
sce.totty.ortho$species <- "Human"
sce.zhuo.ortho$species <- "Rat"

# subset coldata to only clusters and species columns
colData(sce.totty.ortho) <- colData(sce.totty.ortho)[, c("clusters", "species")]
colData(sce.zhuo.ortho) <- colData(sce.zhuo.ortho)[, c("clusters", "species")]

# drop reduced dims
reducedDims(sce.totty.ortho) <- NULL
reducedDims(sce.zhuo.ortho) <- NULL

# subset to human excitatory neurons

sce.combined <- cbind(sce.totty.ortho, sce.zhuo.ortho)
sce.combined
# class: SingleCellExperiment 
# dim: 12069 44568 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(12069): LRP11 PCMT1 ... COL23A1 POMP
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(44568): AAACCCAAGCTAAATG-1 AAACCCACAGGTCCCA-1 ...
#   TTTGTTGGTTGCTCAA-1_19 TTTGTTGTCTTCTGGC-1_19
# colData names(2): clusters species
# reducedDimNames(0):
# mainExpName: RNA
# altExpNames(0):


sce.combined$species <- as.factor(sce.combined$species)
sce.combined$clusters <- as.factor(sce.combined$clusters)

# ====== Copy necessary function from Metaneighbor =======

order_sym_matrix <- function(M, na_value = 0) {
    M <- (M + t(M))/2
    M[is.na(M)] <- na_value
    result <- stats::as.dendrogram(
        stats::hclust(stats::as.dist(1-M), method = "average")
    )
    return(result)
}

# ====== Inhibitory cells ======
#. set up colors
colors <- pals::cols25()[1:length(unique(sce.combined$clusters))]

celltype_colors <- setNames(colors, unique(sce.combined$clusters))

# === Run US Metaneighbor ===

var_genes = variableGenes(dat = sce.combined, exp_labels = sce.combined$species)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = sce.combined,
                             study_id = sce.combined$species,
                             cell_type = sce.combined$clusters,
                             fast_version = TRUE,
                             #one_vs_best=TRUE,
                             #symmetric_output=FALSE
                             )


#pdf(here(plot_dir,"meta_clusters.pdf"))
#plotMetaClusters(mclusters, celltype_NV)
#dev.off()

# === Heatmap Generation ===

# set color scale
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(100))
breaks = seq(0, 1, length=101)
ordering <- order_sym_matrix(celltype_NV)

# --- Heatmap Annotations ---
# get human, baboon, and macaque labels from the row names of celltype_NV
species_anno <- gsub("^(.*?)\\|.*$", "\\1", rownames(celltype_NV))

# get the 3rd rows celltype (CARTPT, PPP1R1B, SST) after | from the row names of celltype_NV
celltype_anno <- gsub("^.*?\\|(.*?)$", "\\1", rownames(celltype_NV))

# get distinct celltype annotation colors
# n <- length(celltype_anno)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]

# make a named vector from celltype_anno and col_vector
#names(col_vector) <- celltype_anno

library(ComplexHeatmap)

# make top (column) annotations
row_ha <- rowAnnotation(
    species = species_anno,
    col = list(
        species = c("Rat" ="#fe9380", "Human" = "#b0d5f5")),
    show_annotation_name = FALSE,
    show_legend = c(TRUE, FALSE)
)

# make row annotations 
top_ha <- HeatmapAnnotation(
        species = species_anno,
    col = list(
        species = c("Rat" ="#fe9380", "Human" = "#b0d5f5")),
    show_annotation_name = FALSE,
    show_legend = c(TRUE, FALSE)
)

# get row order from heatmap to rename labels
hm <- ComplexHeatmap::Heatmap(celltype_NV,
                              cluster_rows = ordering,
                              cluster_columns = ordering,
                              show_column_dend = FALSE, 
                              show_row_dend = FALSE,
                              right_annotation=row_ha,
                              show_column_names=FALSE,
                              row_names_gp = grid::gpar(fontsize = 8),
                              name="AUROC",
                              show_row_names=FALSE
)
row_name_order <- row_order(hm)

# reorder celltype_anno based on row name order
celltype_labels <- celltype_anno[row_name_order]

# # new row labels (only 1 for each celltype in the middle row)
# new_row_labels <- rep("", length(celltype_labels))
# for (i in seq(2, length(celltype_labels), 3)) {
#     new_row_labels[i] <- celltype_labels[i]
# }

# # revert back to original row name order
# new_row_labels <- new_row_labels[order(row_name_order)]

# plot complex heatmap
pdf(here(plot_dir, "Metaneighbor_heatmap_Totty_vs_Zhuo.pdf"), width=6.5, height=5)
ComplexHeatmap::Heatmap(celltype_NV,
                        cluster_rows = ordering,
                        cluster_columns = ordering,
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        right_annotation=row_ha,
                        #top_annotation=top_ha,
                        show_column_names=FALSE,
                        row_names_gp = grid::gpar(fontsize = 8),
                        name="AUROC",
                        #row_labels = new_row_labels,
                        column_title = "Cross-species Mapping of Excitatory Neurons",
                        column_title_side="top",
                        row_title=""
)
dev.off()





# ======== Correlation using spatialLIBD ========

unique(sce.zhuo.ortho$lieden_k25)
# [1] Ppp1r1b_2 Excit_2   Rspo2_2   Rspo2_1   Rspo2_3   Excit_3   Ppp1r1b_3
# [8] Ppp1r1b_1 Excit_1  
# 9 Levels: Ppp1r1b_1 Ppp1r1b_2 Ppp1r1b_3 Rspo2_1 Rspo2_2 Rspo2_3 ... Excit_3

# rename rat lieden clusters to Excit.1 - Excit.9 - order does not matter
# Create a named vector for the new cluster names
new_names <- setNames(paste0("Excit.", seq_along(unique(sce.zhuo.ortho$lieden_k25))),
                      unique(sce.zhuo.ortho$lieden_k25))

# Rename the clusters
sce.zhuo.ortho$lieden_k25 <- as.character(sce.zhuo.ortho$lieden_k25)
sce.zhuo.ortho$lieden_k25 <- new_names[sce.zhuo.ortho$lieden_k25]
sce.zhuo.ortho$lieden_k25 <- factor(sce.zhuo.ortho$lieden_k25, levels = paste0("Excit.", 1:9))

# check names
unique(sce.zhuo.ortho$lieden_k25)
# [1] Excit.1 Excit.2 Excit.3 Excit.4 Excit.5 Excit.6 Excit.7 Excit.8 Excit.9
# 9 Levels: Excit.1 Excit.2 Excit.3 Excit.4 Excit.5 Excit.6 Excit.7 ... Excit.9



colnames(colData(sce.totty.ortho))
#  [1] "orig.ident"             "nCount_originalexp"     "nFeature_originalexp"  
#  [4] "Sample_num"             "Sample"                 "Species"               
#  [7] "Subject"                "Sex"                    "Region"                
# [10] "Subregion"              "DV_axis"                "PI.NeuN"               
# [13] "batch"                  "Barcode"                "sum"                   
# [16] "detected"               "subsets_Mito_sum"       "subsets_Mito_detected" 
# [19] "subsets_Mito_percent"   "total"                  "high_mito"             
# [22] "low_lib"                "low_genes"              "discard_auto"          
# [25] "doubletScore"           "sizeFactor"             "species"               
# [28] "unintegrated_clusters"  "seurat_clusters"        "integrated_snn_res.0.5"
# [31] "ident"                  "key"                    "broad_celltype"        
# [34] "fine_celltype"      

colnames(colData(sce.zhuo.ortho))
#  [1] "orig.ident"                "nCount_RNA"               
#  [3] "nFeature_RNA"              "sample"                   
#  [5] "treatment"                 "addiction.index"          
#  [7] "label"                     "percent.mt"               
#  [9] "nCount_SCT"                "nFeature_SCT"             
# [11] "rfid"                      "integrated_snn_res.0.8"   
# [13] "seurat_clusters"           "cocaine.low"              
# [15] "cocaine.high"              "X933000320047328"         
# [17] "X933000120138592"          "X933000120138586"         
# [19] "X933000320046084"          "X933000320046077"         
# [21] "X933000120138609"          "X933000320186802"         
# [23] "X933000320047225"          "X933000320046609"         
# [25] "X933000320047001"          "X933000320047132"         
# [27] "X933000320186801"          "X933000320046621"         
# [29] "A_933000320046625_JB_257"  "X933000320047104"         
# [31] "X933000320045674"          "Rat_Opioid_HS_1"          
# [33] "Rat_Opioid_HS_2"           "Rat_Amygdala_787A_all_seq"
# [35] "batch"                     "ident"                    
# [37] "sizeFactor"                "lieden_k20"               
# [39] "lieden_k25"                "lieden_k30" 

human_modeling_results <- registration_wrapper(sce.totty.ortho, 
                                                var_registration = "fine_celltype",
                                                var_sample_id = "Sample",
                                                gene_ensembl= "gene_name",
                                                gene_name="gene_name")

rat_modeling_results <- registration_wrapper(sce.zhuo.ortho,
                                                var_registration = "lieden_k25",
                                                var_sample_id = "sample",
                                                gene_ensembl= "ortholog",
                                                gene_name="ortholog")

## extract t-statics and rename
human_t_stats <- human_modeling_results$enrichment[, grep("^t_stat", colnames(human_modeling_results$enrichment))]
colnames(human_t_stats) <- gsub("^t_stat_", "", colnames(human_t_stats))

rat_t_stats <- rat_modeling_results$enrichment[, grep("^t_stat", colnames(rat_modeling_results$enrichment))]
colnames(rat_t_stats) <- gsub("^t_stat_", "", colnames(rat_t_stats))


cor_layer <- layer_stat_cor(
    stats = rat_t_stats,
    modeling_results = human_modeling_results,
    model_type = "enrichment",
    top_n = 500
)

head(cor_layer)
#         ADARB2_TRPS1   ESR1_ADRA1A   GRIK3_TNS3  GULP1_TRHDE MEIS1_PARD3B
# Excit.2   0.16589426  0.1277286277 -0.096328838  0.230882847  -0.05266454
# Excit.8   0.17081398  0.2115518009 -0.051333962 -0.049900702  -0.08631096
# Excit.7   0.09652396 -0.1977296359 -0.019550963 -0.028868341   0.18014852
# Excit.1  -0.14706259 -0.1317785555 -0.040139099 -0.037661865   0.14616610
# Excit.6  -0.15387464  0.0055979794  0.005009133 -0.018260276  -0.02883920
# Excit.3  -0.08231308  0.0022960462  0.043493025  0.051942278  -0.16414986
# Excit.4  -0.00587589  0.0506783209  0.067681953 -0.129260731  -0.05544397
# Excit.5   0.03141996 -0.0007040789  0.055911156 -0.002298512   0.06151763



library(pheatmap)
data_matrix <- cor_layer

# Find the minimum and maximum of your data
min_val <- min(data_matrix)
max_val <- max(data_matrix)

lim <- max(abs(c(min_val, max_val)))

# Because you're scaling rows, it's better to set your breaks based on a standard scale
# like the range from -1 to 1 if your data was normalized or other appropriate values.
breaks_val <- seq(-lim, lim, length.out = 101)  # one more break than the number of colors
my.col <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdBu")))(length(breaks_val))



# Now, you create the heatmap and save it as a PDF
pdf(here(plot_dir, "KeriR01_SpatialRegistration.pdf"), width = 5, height = 5)

pheatmap(t(data_matrix),
         color = my.col,
         #breaks = breaks_val,  # setting breaks to match color transitions
         show_rownames = T,
         show_colnames = T,
         #cutree_cols = 3,
         fontsize=13,
         treeheight_row = 0,
         treeheight_col= 0
         #annotation_col = column_annotations, # Add column annotations
         #annotation_colors = annotation_colors # Specify colors for annotations
)
dev.off()
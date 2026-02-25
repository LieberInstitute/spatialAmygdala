suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("scater")
    library("ggplot2")
    library("patchwork")
})


# load
sce <- readRDS(here("processed-data/snRNAseq/01_build_sce/WHB_amygdala_ITC_subset_SCE.rds"))
sce
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


# ====== Converting from Ensembl to gene symbols ======
# check rowData
head(rowData(sce))
# DataFrame with 6 rows and 1 column
#                 gene_symbol
#                    <factor>
# ENSG00000000003    TSPAN6  
# ENSG00000000005    TNMD    
# ENSG00000000419    DPM1    
# ENSG00000000457    SCYL3   
# ENSG00000000460    C1orf112
# ENSG00000000938    FGR 

# mke rownames gene_id
rowData(sce)$gene_id <- rownames(sce)

# switch to gene_symbol as rownames
rownames(sce) <- rowData(sce)$gene_symbol


# ====== Find sample / library ID for batch correction ======

colnames(colData(sce))
#  [1] "cell_barcode"               "library_label"             
#  [3] "anatomical_division_label"  "index"                     
#  [5] "cell_label"                 "cell_barcode"              
#  [7] "barcoded_cell_sample_label" "library_label"             
#  [9] "feature_matrix_label"       "entity"                    
# [11] "brain_section_label"        "library_method"            
# [13] "donor_label"                "donor_sex"                 
# [15] "dataset_label"              "x"                         
# [17] "y"                          "cluster_alias"             
# [19] "region_of_interest_label"   "anatomical_division_label" 
# [21] "abc_sample_id"              "subcluster"                
# [23] "cluster"                    "supercluster"              
# [25] "neurotransmitter"           "subcluster_color"          
# [27] "cluster_color"              "supercluster_color"        
# [29] "neurotransmitter_color"     "region_of_interest_color" 


length(unique(sce$library_label))
# [1] 572

length(unique(sce$donor_label))
# [1] 4

unique(sce$anatomical_division_label)

# get number of cells per anatomical division 
table(sce$anatomical_division_label)
# Amygdaloid complex    Basal forebrain       Basal nuclei         Cerebellum 
#              15202               6568               7825               1149 
#    Cerebral cortex          Claustrum  Extended amygdala        Hippocampus 
#              16488                833               4378               2634 
#       Hypothalamus           Midbrain     Myelencephalon               Pons 
#               4438               1329                540               1144 
#        Spinal cord           Thalamus 
#                 53               3491 


table(sce$supercluster)
    #       Amygdala excitatory         Cerebellar inhibitory 
    #                       125                           297 
    #           CGE interneuron Eccentric medium spiny neuron 
    #                      4577                         38726 
    # LAMP5-LHX6 and Chandelier             Lower rhombic lip 
    #                       908                          1424 
    #           Mammillary body           Medium spiny neuron 
    #                       389                           117 
    #           MGE interneuron   Midbrain-derived inhibitory 
    #                      4108                          2605 
    #             Miscellaneous                      Splatter 
    #                         3                          6996 
    #       Thalamic excitatory             Upper rhombic lip 
    #                      1825                          3972 

# get number of cells per library_label
table(sce$donor_label)
# H18.30.001 H18.30.002 H19.30.001 H19.30.002 
#        297      23143      21489      21143 


# ======= Normalize, get HVGs, PCA, UMAP ==========

set.seed(1234)
sce <- scater::logNormCounts(sce, assay.type = "counts", log = TRUE, pseudo.count = 1)

# identify highly variable genes
dec <- scran::modelGeneVar(sce, assay.type = "logcounts")
hvg <- scran::getTopHVGs(dec, n = 2000)

# PCA
sce <- scater::runPCA(sce, subset_row = hvg)

# UMAP
sce <- scater::runUMAP(sce, dimred = "PCA", min_dist = 0.3)



# ======= Plot UMAPs of Donor, Anatomical Division, and supercluster =======

pdf(here("plots/snRNAseq/BICCN/WHB_amygdala_ITC_subset_UMAPs.pdf"), width = 18, height = 6)
p1 <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "donor_label", point_size=0.5) +
    ggtitle("UMAP colored by Donor") +
    theme_minimal()

p2 <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "anatomical_division_label", point_size=0.5) +
    ggtitle("UMAP colored by Anatomical Division") +
    theme_minimal()

p3 <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "supercluster", point_size=0.5) +
    ggtitle("UMAP colored by Supercluster") +
    theme_minimal()

p1+p2+p3
dev.off()

# plot ITC markers, TSHZ1, FOXP2, DRD1, GAD1, CPNE4, PRKG1
pdf(here("plots/snRNAseq/BICCN/WHB_amygdala_ITC_subset_UMAPs_ITC_markers.pdf"), width = 18, height = 18)
markers <- c("TSHZ1", "FOXP2", "DRD1", "GAD1", "CPNE4", "PRKG1", "CHST9", "GRM8", "OPRM1")
plot_list <- list()
for (gene in markers) {
    p <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = gene, point_size=0.5) +
        ggtitle(paste("UMAP colored by", gene)) +
        theme_minimal()
    plot_list[[gene]] <- p
}

wrap_plots(plot_list, ncol = 3)
dev.off()

# plot by cluster and  subcluster
pdf(here("plots/snRNAseq/BICCN/WHB_amygdala_ITC_subset_UMAPs_clusters.pdf"), width = 18, height = 6)
p4 <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "cluster", point_size=0.5) +
    ggtitle("UMAP colored by Cluster") +
    theme_minimal()

p4 
dev.off()



# ======== Integrate across donors using Harmony =========

suppressPackageStartupMessages({
    library("harmony")
})

set.seed(1234)
sce_harmony <- sce
harmony_embeddings <- HarmonyMatrix(data = reducedDim(sce_harmony, "PCA"),
                                     meta_data = as.data.frame(colData(sce_harmony)),
                                     vars_use = "donor_label",
                                     do_pca = FALSE)

# store harmony embeddings in reducedDims
reducedDim(sce_harmony, "harmony") <- harmony_embeddings

# UMAP on harmony embeddings
set.seed(1234)
sce_harmony <- scater::runUMAP(sce_harmony, dimred = "harmony", min_dist = 0.3, name = "UMAP_harmony")



# ======= Plot UMAPs of Donor, Anatomical Division, and supercluster after Harmony =======

pdf(here("plots/snRNAseq/BICCN/WHB_amygdala_ITC_subset_UMAPs_harmony.pdf"), width = 18, height = 6)
p1_harmony <- scater::plotReducedDim(sce_harmony, dimred = "UMAP_harmony", colour_by = "donor_label", point_size=0.5) +
    ggtitle("UMAP (Harmony) colored by Donor") +
    theme_minimal() 

p2_harmony <- scater::plotReducedDim(sce_harmony, dimred = "UMAP_harmony", colour_by = "anatomical_division_label", point_size=0.5) +
    ggtitle("UMAP (Harmony) colored by Anatomical Division") +
    theme_minimal()

p3_harmony <- scater::plotReducedDim(sce_harmony, dimred = "UMAP_harmony", colour_by = "supercluster", point_size=0.5) +
    ggtitle("UMAP (Harmony) colored by Supercluster") +
    theme_minimal()

p1_harmony+p2_harmony+p3_harmony
dev.off()


# plot ITC markers, TSHZ1, FOXP2, DRD1, GAD1, CPNE4, PRKG1
pdf(here("plots/snRNAseq/BICCN/WHB_amygdala_ITC_subset_UMAPs_harmony_ITC_markers.pdf"), width = 18, height = 18)
markers <- c("TSHZ1", "FOXP2", "DRD1", "GAD1", "CPNE4", "PRKG1", "CHST9", "GRM8", "OPRM1")
plot_list_harmony <- list()
for (gene in markers) {
    p_harmony <- scater::plotReducedDim(sce_harmony, dimred =   "UMAP_harmony", colour_by = gene, point_size=0.5) +
        ggtitle(paste("UMAP (Harmony) colored by", gene)) +
        theme_minimal()
    plot_list_harmony[[gene]] <- p_harmony
}

wrap_plots(plot_list_harmony, ncol = 3)
dev.off()



# ========== Leiden Clustering on Harmony PCA ==========

suppressPackageStartupMessages({
    library("igraph")
    library("bluster")
    library("scran")
})

# Louvain clustering on harmony-corrected PCA
sce_harmony$louvain_cluster <- clusterCells(
    sce_harmony, 
    use.dimred = "harmony",
    BLUSPARAM = NNGraphParam(
        k = 40,
        cluster.fun = "louvain"
    )
)

pdf(here("plots/snRNAseq/BICCN/WHB_amygdala_ITC_subset_UMAPs_harmony_louvain_clusters.pdf"), width = 6, height = 6)
p_louvain <- scater::plotReducedDim(sce_harmony, dimred = "UMAP_harmony", colour_by = "louvain_cluster", point_size=0.5) +
    ggtitle("UMAP (Harmony) colored by Louvain Clusters") +
    theme_minimal()

p_louvain
dev.off()

# plot but show cluster ID on UMAP
library(ggrepel)
pdf(here("plots/snRNAseq/BICCN/WHB_amygdala_ITC_subset_UMAPs_harmony_louvain_clusters_labels.pdf"), width = 6, height = 6)
umap_df <- as.data.frame(reducedDim(sce_harmony, "UMAP_harmony"))
umap_df$louvain_cluster <- sce_harmony$louvain_cluster  

p_louvain_labels <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = louvain_cluster)) +
    geom_point(size = 0.5, alpha = 0.6) +
    theme_minimal() +
    ggtitle("UMAP (Harmony) colored by Louvain Clusters") +
    geom_text_repel(data = aggregate(cbind(UMAP1, UMAP2) ~ louvain_cluster, umap_df, mean),
                    aes(label = louvain_cluster),
                    size = 5,
                    color = "black")
p_louvain_labels
dev.off()



# ========= Subset True ITCS ==========

true_itcs <- c("8", "10", "15", "16", "17", "18", "19")

# subset
sce_itcs <- sce_harmony[, sce_harmony$louvain_cluster %in% true_itcs]
sce_itcs
# class: SingleCellExperiment 
# dim: 59357 39276 
# metadata(0):
# assays(2): counts logcounts
# rownames(59357): TSPAN6 TNMD ... AC114982.3 AC084756.2
# rowData names(2): gene_symbol gene_id
# colnames(39276): 10X358_7:TATACCTTCATCGCAA 10X190_1:GACCAATGTTTCGGCG
#   ... 10X349_4:GACTCAACAGACCGCT 10X354_7:CGCATGGGTCAGTCGC
# colData names(32): cell_barcode library_label ... sizeFactor
#   louvain_cluster
# reducedDimNames(4): PCA UMAP harmony UMAP_harmony
# mainExpName: NULL
# altExpNames(0):

# ========== Redo normalization, HVGs, PCA, UMAP on ITCs only ==========

# normalize
sce_itcs <- scater::logNormCounts(sce_itcs, assay.type = "counts", log = TRUE, pseudo.count = 1)

# identify highly variable genes
dec_itcs <- scran::modelGeneVar(sce_itcs, assay.type = "logcounts")
hvg_itcs <- scran::getTopHVGs(dec_itcs, n = 2000)

# PCA
sce_itcs <- scater::runPCA(sce_itcs, subset_row = hvg_itcs)

# Harmony

sce_itcs_harmony <- sce_itcs
harmony_embeddings_itcs <- HarmonyMatrix(data = reducedDim(sce_itcs_harmony, "PCA"),
                                     meta_data = as.data.frame(colData(sce_itcs_harmony)),
                                     vars_use = "donor_label",
                                     do_pca = FALSE)

# store harmony embeddings in reducedDims
reducedDim(sce_itcs_harmony, "harmony") <- harmony_embeddings_itcs

# UMAP on harmony embeddings
sce_itcs_harmony <- scater::runUMAP(sce_itcs_harmony, dimred = "harmony", min_dist = 0.3, name = "UMAP_harmony")


# ======= Plot UMAPs of Donor, Anatomical Division, and supercluster after Harmony on ITCs =======

pdf(here("plots/snRNAseq/BICCN/WHB_amygdala_ITC_subset_ITCs_UMAPs_harmony.pdf"), width = 18, height = 6)
p1_itcs_harmony <- scater::plotReducedDim(sce_itcs_harmony, dimred = "UMAP_harmony", colour_by = "donor_label", point_size=0.5) +
    ggtitle("UMAP (Harmony) colored by Donor") +
    theme_minimal()

p2_itcs_harmony <- scater::plotReducedDim(sce_itcs_harmony, dimred = "UMAP_harmony", colour_by = "anatomical_division_label", point_size=0.5) +
    ggtitle("UMAP (Harmony) colored by Anatomical Division") +
    theme_minimal()

p3_itcs_harmony <- scater::plotReducedDim(sce_itcs_harmony, dimred = "UMAP_harmony", colour_by = "supercluster", point_size=0.5) +
    ggtitle("UMAP (Harmony) colored by Supercluster") +
    theme_minimal()

p1_itcs_harmony+p2_itcs_harmony+p3_itcs_harmony
dev.off()

# plot ITC markers, TSHZ1, FOXP2, DRD1, GAD1, CPNE4, PRKG1
pdf(here("plots/snRNAseq/BICCN/WHB_amygdala_ITC_subset_ITCs_UMAPs_harmony_ITC_markers.pdf"), width = 18, height = 18)
markers <- c("TSHZ1", "FOXP2", "DRD1", "GAD1", "CPNE4", "PRKG1", "CHST9", "GRM8", "OPRM1")
plot_list_itcs_harmony <- list()
for (gene in markers) {
    p_itcs_harmony <- scater::plotReducedDim(sce_itcs_harmony, dimred =   "UMAP_harmony", colour_by = gene, point_size=0.5) +
        ggtitle(paste("UMAP (Harmony) colored by", gene)) +
        theme_minimal()
    plot_list_itcs_harmony[[gene]] <- p_itcs_harmony
}

wrap_plots(plot_list_itcs_harmony, ncol = 3)
dev.off()


# cluster
sce_itcs_harmony$louvain_cluster <- clusterCells(
    sce_itcs_harmony, 
    use.dimred = "harmony",
    BLUSPARAM = NNGraphParam(
        k = 50,
        cluster.fun = "louvain"
    )
)

# plot clusters
pdf(here("plots/snRNAseq/BICCN/WHB_amygdala_ITC_subset_ITCs_UMAPs_harmony_louvain_clusters.pdf"), width = 6, height = 6)
p_louvain_itcs <- scater::plotReducedDim(sce_itcs_harmony, dimred = "UMAP_harmony", colour_by = "louvain_cluster", point_size=0.5) +
    ggtitle("UMAP (Harmony) colored by Louvain Clusters") +
    theme_minimal() 
p_louvain_itcs
dev.off()


# find one vs all markers and save markers table
# Option 2: Specify direction in findMarkers (re-run)
markers_itcs <- findMarkers(
    sce_itcs_harmony, 
    groups = sce_itcs_harmony$louvain_cluster, 
    assay.type = "logcounts",
    direction = "up"
)

library(dplyr)
library(purrr)

# Option 1: Use bind_rows (fills missing with NA)
top_markers_df <- imap_dfr(markers_itcs, function(df, cluster_name) {
    df |>
        as.data.frame() |>
        head(100) |>
        mutate(
            gene = rownames(df)[1:min(100, nrow(df))],
            cluster = cluster_name
        )
}) |>
    relocate(cluster, gene)

write.csv(top_markers_df, here("processed-data/snRNAseq/BICCN/top100_markers_per_cluster.csv"), row.names = FALSE)

# save final ITC SCE
saveRDS(sce_itcs_harmony, here("processed-data/snRNAseq/BICCN/WHB_amygdala_Final_ITC_subsets.rds"))
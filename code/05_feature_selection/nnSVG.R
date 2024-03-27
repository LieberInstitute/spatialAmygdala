# qrsh -pe local 8 -l mem_free=15G,h_vmem=15G -now n

library("SpatialExperiment")
library("here")
library("ggspavis")
library("scater")
library("pheatmap")
library("spatialLIBD")
library("patchwork")
library("nnSVG")
library("scran")
library("BayesSpace")

# save directiories
plot_dir = here("plots", "05_feature_selection")
processed_dir = here("processed-data", "05_feature_selection")

# load object
load(here("processed-data","03_qc_metrics","spe_local_outliers.Rdata"))
spe
# class: SpatialExperiment 
# dim: 28412 29885 
# metadata(0):
#     assays(1): counts
# rownames(28412): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(6): source type ... gene_name Symbol.uniq
# colnames(29885): AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 ...
# TTGTTTCATTAGTCTA-1 TTGTTTCCATACAACT-1
# colData names(42): sample_id in_tissue ... qc_detected discard
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
#     spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


# -------- Spatially aware feature selection -------
# because we have multiple samples (visium arrays) we will need to run nnSVG
# once per sample and then average the top SVGs, according the to the vingette

# check sample IDs
table(colData(spe)$sample_id)
# V13M06-387_A1 V13M06-387_B1 V13M06-387_C1 V13M06-387_D1 V13M06-388_A1 
# 3694          3882          3453          3280          3986 
# V13M06-388_B1 V13M06-388_C1 V13M06-388_D1 
# 3571          3999          4020 


# run nnSVG once per sample and store lists of top SVGs

sample_ids <-unique(colData(spe)$sample_id)

res_list <- as.list(rep(NA, length(sample_ids)))
names(res_list) <- sample_ids

for (s in seq_along(sample_ids)) {
    
    Sys.time()
    # select sample
    ix <- colData(spe)$sample_id == sample_ids[s]
    spe_sub <- spe[, ix]
    
    dim(spe_sub)
    
    # run nnSVG filtering for mitochondrial genes and low-expressed genes
    # (note: set 'filter_mito = TRUE' in most datasets)
    spe_sub <- filter_genes(spe_sub)
    
    # remove any zeros introduced by filtering
    ix_zeros <- colSums(counts(spe_sub)) == 0
    if (sum(ix_zeros) > 0) {
        spe_sub <- spe_sub[, !ix_zeros]
    }
    
    dim(spe_sub)
    
    # re-calculate logcounts after filtering
    spe_sub <- computeLibraryFactors(spe_sub)
    spe_sub <- logNormCounts(spe_sub)
    
    # run nnSVG
    set.seed(123)
    spe_sub <- nnSVG(spe_sub)
    
    # store results for this sample
    res_list[[s]] <- rowData(spe_sub)
}

# number of genes that passed filtering (and subsampling) for each sample
sapply(res_list, nrow)
# V13M06-387_A1 V13M06-387_B1 V13M06-387_C1 V13M06-387_D1 V13M06-388_A1 
# 30            25            13            18            28 
# V13M06-388_B1 V13M06-388_C1 V13M06-388_D1 
# 30            25            27 

# match results from each sample and store in matching rows
res_ranks <- matrix(NA, nrow = nrow(spe), ncol = length(sample_ids))
rownames(res_ranks) <- rownames(spe)
colnames(res_ranks) <- sample_ids

for (s in seq_along(sample_ids)) {
    stopifnot(colnames(res_ranks)[s] == sample_ids[s])
    stopifnot(colnames(res_ranks)[s] == names(res_list)[s])
    
    rownames_s <- rownames(res_list[[s]])
    res_ranks[rownames_s, s] <- res_list[[s]][, "rank"]
}


# remove genes that were filtered out in all samples
ix_allna <- apply(res_ranks, 1, function(r) all(is.na(r)))
res_ranks <- res_ranks[!ix_allna, ]

dim(res_ranks)
# [1] 30  8


# calculate average ranks
# note missing values due to filtering for samples
avg_ranks <- rowMeans(res_ranks, na.rm = TRUE)


# calculate number of samples where each gene is within top 100 ranked SVGs
# for that sample
n_withinTop100 <- apply(res_ranks, 1, function(r) sum(r <= 100, na.rm = TRUE))

# summary table
df_summary <- data.frame(
    gene_id = names(avg_ranks), 
    gene_name = rowData(spe)[names(avg_ranks), "gene_name"], 
    gene_type = rowData(spe)[names(avg_ranks), "type"], 
    overall_rank = rank(avg_ranks), 
    average_rank = unname(avg_ranks), 
    n_withinTop100 = unname(n_withinTop100), 
    row.names = names(avg_ranks)
)

# sort by average rank
df_summary <- df_summary[order(df_summary$average_rank), ]
write.csv(df_summary, here(processed_dir, "nnSVG_summary.csv"), row.names=FALSE)
head(df_summary)
# gene_id gene_name gene_type overall_rank average_rank n_withinTop100
# SH3GL2  SH3GL2    SH3GL2      gene            1     2.125000              8
# CD44      CD44      CD44      gene            2     3.666667              3
# TAC1      TAC1      TAC1      gene            3     4.250000              8
# RAB3C    RAB3C     RAB3C      gene            4     5.000000              8
# RPRM      RPRM      RPRM      gene            5     5.666667              3
# ABR        ABR       ABR      gene            6     6.000000              8

ix_gene <- which(rowData(spe)$gene_name == df_summary[1, "gene_name"])

df <- as.data.frame(cbind(
    colData(spe), 
    spatialCoords(spe), 
    gene = counts(spe)[ix_gene, ]
))

p<-ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = gene)) + 
    facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
    geom_point(size = 0.6) + 
    scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                         name = "counts", breaks = range(df$gene)) + 
    scale_y_reverse() + 
    ggtitle(paste0(rowData(spe)$gene_name[ix_gene], " expression")) + 
    theme_bw() + 
    theme(aspect.ratio = 1, 
          panel.grid = element_blank(), 
          plot.title = element_text(face = "italic"), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())

pdf(here(plot_dir,"TopSVG_spotplot.pdf"))
print(p)
dev.off()
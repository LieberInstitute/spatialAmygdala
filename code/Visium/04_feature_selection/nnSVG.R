library("SpatialExperiment")
library("here")
library("ggspavis")
library("scater")
library("spatialLIBD")
library("patchwork")
library("nnSVG")
library("scran")

# save directiories
plot_dir = here("plots", "05_feature_selection")
processed_dir = here("processed-data", "05_feature_selection")

# load object
load(here("processed-data","04_normalization","spe_stitched_norm.Rdata"))
spe
# class: SpatialExperiment 
# dim: 36601 227302 
# metadata(0):
# assays(2): counts logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
#   ENSG00000277196
# rowData names(1): symbol
# colnames(227302): AAACAAGTATCTCCCA-1_V13Y24-346_A1
#   AAACAATCTACTAGCA-1_V13Y24-346_A1 ... TTGTTTCATTAGTCTA-1_V13B23-407_C1
#   TTGTTTCCATACAACT-1_V13B23-407_C1
# colData names(37): in_tissue array_row ... subsets_mito_percent_z
#   sizeFactor
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


# -------- Spatially aware feature selection -------
# because we have multiple samples (arrays combined into donors) we will need to run nnSVG
# once per sample and then average the top SVGs, according the to the vingette

# check sample IDs
table(colData(spe)$sample_id)
# Br2743 Br6423 Br6471 Br6660 Br8325 Br9017 Br9192 Br9206 Br9280 Br9469 
#  32855  34353  37795  37085  29885  26382  28608  26850  26722  36639 


# run nnSVG once per sample and store lists of top SVGs
sample_ids <-unique(colData(spe)$sample_id)

rownames(spe) <- rowData(spe)$symbol
rowData(spe)$gene_name <- rowData(spe)$symbol

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
    gene_name = rowData(spe)[names(avg_ranks), "symbol"], 
    overall_rank = rank(avg_ranks), 
    average_rank = unname(avg_ranks), 
    n_withinTop100 = unname(n_withinTop100), 
    row.names = names(avg_ranks)
)

# sort by average rank
df_summary <- df_summary[order(df_summary$average_rank), ]
write.csv(df_summary, here(processed_dir, "nnSVG_summary_stitched.csv"), row.names=FALSE)
head(df_summary)
# gene_id gene_name gene_type overall_rank average_rank n_withinTop100
# SH3GL2  SH3GL2    SH3GL2      gene            1     2.125000              8
# CD44      CD44      CD44      gene            2     3.666667              3
# TAC1      TAC1      TAC1      gene            3     4.250000              8
# RAB3C    RAB3C     RAB3C      gene            4     5.000000              8
# RPRM      RPRM      RPRM      gene            5     5.666667              3
# ABR        ABR       ABR      gene            6     6.000000              8

ix_gene <- which(rowData(spe)$symbol == df_summary[1, "gene_name"])

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
    ggtitle(paste0(rowData(spe)$symbol[ix_gene], " expression")) + 
    theme_bw() + 
    theme(aspect.ratio = 1, 
          panel.grid = element_blank(), 
          plot.title = element_text(face = "italic"), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())

pdf(here(plot_dir,"TopSVG_spotplot_stitched.pdf"))
print(p)
dev.off()
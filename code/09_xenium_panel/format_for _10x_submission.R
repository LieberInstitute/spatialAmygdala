suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("Seurat")
})

# Save directories
processed_dir = here("processed-data","09_xenium_panel")

load(here("processed-data","08_clustering", "spe_clusters_k10.Rdata"), verbose = TRUE)
spe
# class: SpatialExperiment 
# dim: 28412 29885 
# metadata(2): BayesSpace.data chain.h5
# assays(2): counts logcounts
# rownames(28412): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(7): source type ... Symbol.uniq is.HVG
# colnames(29885): AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 ... TTGTTTCATTAGTCTA-1 TTGTTTCCATACAACT-1
# colData names(47): sample_id in_tissue ... cluster.init spatial.cluster
# reducedDimNames(4): 10x_pca 10x_tsne 10x_umap PCA
# mainExpName: NULL
# altExpNames(0):
#     spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

# drop duplicate colData
colData(spe) <- colData(spe)[ , !duplicated(colnames(colData(spe)))]

# add names to spatial clusters
spe$spatial.cluster <- factor(spe$spatial.cluster)
levels(spe$spatial.cluster) <- c('WM.1', 'aBA', 'LA.1', 'BA', 'vmBA', 'GABA.n', 'LA.2', 'EC', 'WM.2', 'Endo')
unique(spe$spatial.cluster)

# === Get required counts, gene_ids, gene_symbol, and barcodes
library("DropletUtils")

dim(spe)
# [1] 28412 29885

colnames(rowData(spe))
# [1] "source"       "type"         "gene_id"      "gene_version" "gene_name"    "Symbol.uniq"  "is.HVG"    

colnames(colData(spe))
# [1] "sample_id"              "in_tissue"              "array_row"              "array_col"             
# [5] "10x_graphclust"         "10x_kmeans_10_clusters" "10x_kmeans_2_clusters"  "10x_kmeans_3_clusters" 
# [9] "10x_kmeans_4_clusters"  "10x_kmeans_5_clusters"  "10x_kmeans_6_clusters"  "10x_kmeans_7_clusters" 
# [13] "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                    "sum_umi"               
# [17] "sum_gene"               "expr_chrM"              "expr_chrM_ratio"        "ManualAnnotation"      
# [21] "slide"                  "array"                  "brnum"                  "species"               
# [25] "replicate"              "overlaps_tissue"        "sum"                    "detected"              
# [29] "subsets_mito_sum"       "subsets_mito_detected"  "subsets_mito_percent"   "total"                 
# [33] "qc_lib_size"            "qc_mito"                "qc_detected"            "discard"               
# [37] "sizeFactor"             "row"                    "col"                    "cluster.init"          
# [41] "spatial.cluster"     



# create MEX files
write10xCounts(
    here("processed-data","09_xenium_panel","xenium_mex_files"),
    counts(spe),
    gene.id = rowData(spe)$gene_id,
    gene.symbol = rowData(spe)$Symbol.uniq,
    barcodes = colData(spe)$key,
    type = "sparse",
    version = "3"
)

list.files(here("processed-data","09_xenium_panel","xenium_mex_files"))
# [1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"



# ====== Adding spatial domain annotations ======
# bundleOutputs is a function provided by 10x. Note that their code had a typo where
# the "data" argument was not actually called and instead directly used "seurat_obj".
# 
# I have corrected this below.

# Define function
bundleOutputs <- function(out_dir, data, barcodes = colnames(data), cell_type = "cell_type", subset = 1:length(barcodes)) {
    
    if (require("data.table", quietly = TRUE)) {
        data.table::fwrite(
            data.table::data.table(
                barcode = barcodes,
                annotation = unlist(data[[cell_type]])
            )[subset, ],
            file.path(out_dir, "annotations.csv")
        )
    } else {
        write.table(
            data.frame(
                barcode = barcodes,
                annotation = unlist(data[[cell_type]])
            )[subset, ],
            file.path(out_dir, "annotations.csv"),
            sep = ",", row.names = FALSE
        )
    }
    
    bundle <- file.path(out_dir, paste0(basename(out_dir), ".zip"))
    
    utils::zip(
        bundle,
        list.files(out_dir, full.names = TRUE),
        zip = "zip"
    )
    
    if (file.info(bundle)$size / 1e6 > 500) {
        warning("The output file is more than 500 MB and will need to be subset further.")
    }
}

# Run function
bundleOutputs(out_dir = here("processed-data","09_xenium_panel","xenium_mex_files"), 
              data = spe, 
              barcodes = spe$key,
              cell_type = "spatial.cluster")

list.files(here("processed-data","09_xenium_panel","xenium_mex_files"))
# [1] "annotations.csv"      "barcodes.tsv.gz"      "features.tsv.gz"      "matrix.mtx.gz"        "xenium_mex_files.zip"


# read annotations.csv to check
test <- read.csv(here("processed-data","09_xenium_panel","xenium_mex_files","annotations.csv"))
head(test)

# check to see if there are any duplicates in barcode
any(duplicated(test$barcode))

session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-10-12 r85331)
# os       macOS Sonoma 14.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Denver
# date     2024-01-27
# rstudio  2023.09.1+494 Desert Sunflower (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version   date (UTC) lib source
# abind                    1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
# AnnotationDbi            1.64.0    2023-10-24 [1] Bioconductor
# AnnotationHub            3.10.0    2023-10-24 [1] Bioconductor
# assertthat               0.2.1     2019-03-21 [1] CRAN (R 4.3.0)
# attempt                  0.3.1     2020-05-03 [1] CRAN (R 4.3.0)
# BayesSpace             * 1.12.0    2023-10-24 [1] Bioconductor
# beachmat                 2.18.0    2023-10-24 [1] Bioconductor
# beeswarm                 0.4.0     2021-06-01 [1] CRAN (R 4.3.0)
# benchmarkme              1.0.8     2022-06-12 [1] CRAN (R 4.3.0)
# benchmarkmeData          1.0.4     2020-04-23 [1] CRAN (R 4.3.0)
# Biobase                * 2.62.0    2023-10-24 [1] Bioconductor
# BiocFileCache            2.11.1    2023-11-03 [1] Github (Bioconductor/BiocFileCache@7b7c4d3)
# BiocGenerics           * 0.48.1    2023-11-01 [1] Bioconductor
# BiocIO                   1.12.0    2023-10-24 [1] Bioconductor
# BiocManager              1.30.22   2023-08-08 [1] CRAN (R 4.3.0)
# BiocNeighbors            1.20.0    2023-10-24 [1] Bioconductor
# BiocParallel             1.36.0    2023-10-24 [1] Bioconductor
# BiocSingular             1.18.0    2023-10-24 [1] Bioconductor
# BiocVersion              3.18.0    2023-05-11 [1] Bioconductor
# Biostrings               2.70.1    2023-10-25 [1] Bioconductor
# bit                      4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
# bit64                    4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
# bitops                   1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# blob                     1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
# bluster                  1.12.0    2023-10-24 [1] Bioconductor
# bslib                    0.5.1     2023-08-11 [1] CRAN (R 4.3.0)
# cachem                   1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
# cli                      3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# cluster                  2.1.4     2022-08-22 [1] CRAN (R 4.3.1)
# coda                     0.19-4    2020-09-30 [1] CRAN (R 4.3.0)
# codetools                0.2-19    2023-02-01 [1] CRAN (R 4.3.1)
# colorspace               2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# config                   0.3.2     2023-08-30 [1] CRAN (R 4.3.0)
# cowplot                  1.1.1     2020-12-30 [1] CRAN (R 4.3.0)
# crayon                   1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# curl                     5.2.0     2023-12-08 [1] CRAN (R 4.3.1)
# data.table               1.14.8    2023-02-17 [1] CRAN (R 4.3.0)
# DBI                      1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
# dbplyr                   2.4.0     2023-10-26 [1] CRAN (R 4.3.1)
# DelayedArray             0.28.0    2023-10-24 [1] Bioconductor
# DelayedMatrixStats       1.24.0    2023-10-24 [1] Bioconductor
# deldir                   1.0-9     2023-05-17 [1] CRAN (R 4.3.0)
# digest                   0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
# DirichletReg             0.7-1     2021-05-18 [1] CRAN (R 4.3.0)
# doParallel               1.0.17    2022-02-07 [1] CRAN (R 4.3.0)
# dotCall64                1.1-0     2023-10-17 [1] CRAN (R 4.3.1)
# dplyr                    1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# dqrng                    0.3.1     2023-08-30 [1] CRAN (R 4.3.0)
# DropletUtils           * 1.22.0    2023-10-24 [1] Bioconductor
# DT                       0.30      2023-10-05 [1] CRAN (R 4.3.1)
# edgeR                    4.0.1     2023-10-30 [1] Bioconductor
# ellipsis                 0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
# ExperimentHub            2.10.0    2023-10-24 [1] Bioconductor
# fansi                    1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# fastmap                  1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# fields                   15.2      2023-08-17 [1] CRAN (R 4.3.0)
# filelock                 1.0.3     2023-12-11 [1] CRAN (R 4.3.1)
# fitdistrplus             1.1-11    2023-04-25 [1] CRAN (R 4.3.0)
# foreach                  1.5.2     2022-02-02 [1] CRAN (R 4.3.0)
# Formula                  1.2-5     2023-02-24 [1] CRAN (R 4.3.0)
# future                   1.33.0    2023-07-01 [1] CRAN (R 4.3.0)
# future.apply             1.11.0    2023-05-21 [1] CRAN (R 4.3.0)
# generics                 0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.38.1    2023-11-11 [1] Bioconductor
# GenomeInfoDbData         1.2.11    2023-10-27 [1] Bioconductor
# GenomicAlignments        1.38.0    2023-10-24 [1] Bioconductor
# GenomicRanges          * 1.54.1    2023-10-30 [1] Bioconductor
# ggbeeswarm               0.7.2     2023-04-29 [1] CRAN (R 4.3.0)
# ggplot2                * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
# ggrepel                  0.9.4     2023-10-13 [1] CRAN (R 4.3.1)
# ggridges                 0.5.4     2022-09-26 [1] CRAN (R 4.3.0)
# globals                  0.16.2    2022-11-21 [1] CRAN (R 4.3.0)
# glue                     1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
# goftest                  1.2-3     2021-10-07 [1] CRAN (R 4.3.0)
# golem                    0.4.1     2023-06-05 [1] CRAN (R 4.3.0)
# gridExtra              * 2.3       2017-09-09 [1] CRAN (R 4.3.0)
# gtable                   0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# HDF5Array                1.30.0    2023-10-24 [1] Bioconductor
# here                   * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# htmltools                0.5.7     2023-11-03 [1] CRAN (R 4.3.1)
# htmlwidgets              1.6.2     2023-03-17 [1] CRAN (R 4.3.0)
# httpuv                   1.6.12    2023-10-23 [1] CRAN (R 4.3.1)
# httr                     1.4.7     2023-08-15 [1] CRAN (R 4.3.0)
# ica                      1.0-3     2022-07-08 [1] CRAN (R 4.3.0)
# igraph                   1.5.1     2023-08-10 [1] CRAN (R 4.3.0)
# interactiveDisplayBase   1.40.0    2023-10-24 [1] Bioconductor
# IRanges                * 2.36.0    2023-10-24 [1] Bioconductor
# irlba                    2.3.5.1   2022-10-03 [1] CRAN (R 4.3.0)
# iterators                1.0.14    2022-02-05 [1] CRAN (R 4.3.0)
# jquerylib                0.1.4     2021-04-26 [1] CRAN (R 4.3.0)
# jsonlite                 1.8.8     2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST                 1.42.0    2023-10-24 [1] Bioconductor
# KernSmooth               2.23-22   2023-07-10 [1] CRAN (R 4.3.1)
# later                    1.3.1     2023-05-02 [1] CRAN (R 4.3.0)
# lattice                  0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# lazyeval                 0.2.2     2019-03-15 [1] CRAN (R 4.3.0)
# leiden                   0.4.3     2022-09-10 [1] CRAN (R 4.3.0)
# lifecycle                1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# limma                    3.58.0    2023-10-26 [1] Bioconductor
# listenv                  0.9.0     2022-12-16 [1] CRAN (R 4.3.0)
# lmtest                   0.9-40    2022-03-21 [1] CRAN (R 4.3.0)
# locfit                   1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
# magick                   2.8.1     2023-10-22 [1] CRAN (R 4.3.1)
# magrittr                 2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# maps                     3.4.1.1   2023-11-03 [1] CRAN (R 4.3.1)
# MASS                     7.3-60    2023-05-04 [1] CRAN (R 4.3.1)
# Matrix                   1.6-1.1   2023-09-18 [1] CRAN (R 4.3.1)
# MatrixGenerics         * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats            * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# maxLik                   1.5-2     2021-07-26 [1] CRAN (R 4.3.0)
# mclust                   6.0.1     2023-11-15 [1] CRAN (R 4.3.1)
# memoise                  2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
# metapod                  1.10.0    2023-10-24 [1] Bioconductor
# mime                     0.12      2021-09-28 [1] CRAN (R 4.3.0)
# miniUI                   0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
# miscTools                0.6-28    2023-05-03 [1] CRAN (R 4.3.0)
# munsell                  0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# nlme                     3.1-163   2023-08-09 [1] CRAN (R 4.3.1)
# paletteer                1.5.0     2022-10-19 [1] CRAN (R 4.3.0)
# parallelly               1.36.0    2023-05-26 [1] CRAN (R 4.3.0)
# patchwork              * 1.1.3     2023-08-14 [1] CRAN (R 4.3.0)
# pbapply                  1.7-2     2023-06-27 [1] CRAN (R 4.3.0)
# pillar                   1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig                2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# plotly                   4.10.3    2023-10-21 [1] CRAN (R 4.3.1)
# plyr                     1.8.9     2023-10-02 [1] CRAN (R 4.3.1)
# png                      0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
# polyclip                 1.10-6    2023-09-27 [1] CRAN (R 4.3.1)
# progressr                0.14.0    2023-08-10 [1] CRAN (R 4.3.0)
# promises                 1.2.1     2023-08-10 [1] CRAN (R 4.3.0)
# purrr                    1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
# R.methodsS3              1.8.2     2022-06-13 [1] CRAN (R 4.3.0)
# R.oo                     1.25.0    2022-06-12 [1] CRAN (R 4.3.0)
# R.utils                  2.12.2    2022-11-11 [1] CRAN (R 4.3.0)
# R6                       2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# RANN                     2.6.1     2019-01-08 [1] CRAN (R 4.3.0)
# rappdirs                 0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer           * 1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                     1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
# RcppAnnoy                0.0.21    2023-07-02 [1] CRAN (R 4.3.0)
# RCurl                    1.98-1.13 2023-11-02 [1] CRAN (R 4.3.1)
# rematch2                 2.1.2     2020-05-01 [1] CRAN (R 4.3.0)
# reshape2                 1.4.4     2020-04-09 [1] CRAN (R 4.3.0)
# restfulr                 0.0.15    2022-06-16 [1] CRAN (R 4.3.0)
# reticulate               1.34.0    2023-10-12 [1] CRAN (R 4.3.1)
# rhdf5                    2.46.0    2023-10-24 [1] Bioconductor
# rhdf5filters             1.14.0    2023-10-24 [1] Bioconductor
# Rhdf5lib                 1.24.0    2023-10-24 [1] Bioconductor
# rjson                    0.2.21    2022-01-09 [1] CRAN (R 4.3.0)
# rlang                    1.1.2     2023-11-04 [1] CRAN (R 4.3.1)
# ROCR                     1.0-11    2020-05-02 [1] CRAN (R 4.3.0)
# rprojroot                2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
# Rsamtools                2.18.0    2023-10-24 [1] Bioconductor
# RSQLite                  2.3.4     2023-12-08 [1] CRAN (R 4.3.1)
# rstudioapi               0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# rsvd                     1.0.5     2021-04-16 [1] CRAN (R 4.3.0)
# rtracklayer              1.62.0    2023-10-24 [1] Bioconductor
# Rtsne                    0.16      2022-04-17 [1] CRAN (R 4.3.0)
# S4Arrays                 1.2.0     2023-10-24 [1] Bioconductor
# S4Vectors              * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# sandwich                 3.0-2     2022-06-15 [1] CRAN (R 4.3.0)
# sass                     0.4.7     2023-07-15 [1] CRAN (R 4.3.0)
# ScaledMatrix             1.10.0    2023-10-24 [1] Bioconductor
# scales                   1.2.1     2022-08-20 [1] CRAN (R 4.3.0)
# scater                   1.30.0    2023-10-24 [1] Bioconductor
# scattermore              1.2       2023-06-12 [1] CRAN (R 4.3.0)
# scran                    1.30.0    2023-10-24 [1] Bioconductor
# sctransform              0.4.1     2023-10-19 [1] CRAN (R 4.3.1)
# scuttle                  1.12.0    2023-10-24 [1] Bioconductor
# sessioninfo            * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# Seurat                 * 4.4.0     2023-09-28 [1] CRAN (R 4.3.1)
# SeuratObject           * 5.0.0     2023-10-26 [1] CRAN (R 4.3.1)
# shiny                    1.7.5.1   2023-10-14 [1] CRAN (R 4.3.1)
# shinyWidgets             0.8.0     2023-08-30 [1] CRAN (R 4.3.0)
# SingleCellExperiment   * 1.24.0    2023-10-24 [1] Bioconductor
# sp                       2.1-1     2023-10-16 [1] CRAN (R 4.3.1)
# spam                     2.10-0    2023-10-23 [1] CRAN (R 4.3.1)
# SparseArray              1.2.2     2023-11-08 [1] Bioconductor
# sparseMatrixStats        1.14.0    2023-10-24 [1] Bioconductor
# SpatialExperiment      * 1.12.0    2023-10-26 [1] Bioconductor
# spatialLIBD            * 1.13.4    2023-05-25 [1] Bioconductor
# spatstat.data            3.0-3     2023-10-24 [1] CRAN (R 4.3.1)
# spatstat.explore         3.2-5     2023-10-22 [1] CRAN (R 4.3.1)
# spatstat.geom            3.2-7     2023-10-20 [1] CRAN (R 4.3.1)
# spatstat.random          3.2-1     2023-10-21 [1] CRAN (R 4.3.1)
# spatstat.sparse          3.0-3     2023-10-24 [1] CRAN (R 4.3.1)
# spatstat.utils           3.0-4     2023-10-24 [1] CRAN (R 4.3.1)
# statmod                  1.5.0     2023-01-06 [1] CRAN (R 4.3.0)
# stringi                  1.8.3     2023-12-11 [1] CRAN (R 4.3.1)
# stringr                  1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment   * 1.32.0    2023-10-24 [1] Bioconductor
# survival                 3.5-7     2023-08-14 [1] CRAN (R 4.3.1)
# tensor                   1.5       2012-05-05 [1] CRAN (R 4.3.0)
# tibble                   3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                    1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
# tidyselect               1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                     1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# uwot                     0.1.16    2023-06-29 [1] CRAN (R 4.3.0)
# vctrs                    0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
# vipor                    0.4.5     2017-03-22 [1] CRAN (R 4.3.0)
# viridis                  0.6.4     2023-07-22 [1] CRAN (R 4.3.0)
# viridisLite              0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
# withr                    2.5.2     2023-10-30 [1] CRAN (R 4.3.1)
# xgboost                  1.7.5.1   2023-03-30 [1] CRAN (R 4.3.0)
# XML                      3.99-0.15 2023-11-02 [1] CRAN (R 4.3.1)
# xtable                   1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
# XVector                  0.42.0    2023-10-24 [1] Bioconductor
# yaml                     2.3.7     2023-01-23 [1] CRAN (R 4.3.0)
# zlibbioc                 1.48.0    2023-10-24 [1] Bioconductor
# zoo                      1.8-12    2023-04-13 [1] CRAN (R 4.3.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


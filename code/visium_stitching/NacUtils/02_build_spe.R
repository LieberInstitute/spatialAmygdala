args <- commandArgs(trailingOnly = TRUE)
group <- args[1]
group

setwd('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/')
library(here)
library(tidyverse)
library(visiumStitched)
library(sessioninfo)
library(Matrix)
library(SpatialExperiment)
library(HDF5Array)
library(scran)
library(scater)
library(scry)
library(BiocParallel)

info_path = here('processed-data', 'visium_stitching', 'NacUtils', paste0(group,'.csv'))
coords_dir = here('processed-data', 'visium_stitching', 'NacUtils', 'spe_inputs')
spe_dir = here('processed-data', 'visium_stitching', 'NacUtils', 'spe')
plot_dir = here('plots', 'visium_stitching', 'NacUtils')

wm_genes = c("MBP", "GFAP")

num_cores = as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
dir.create(coords_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Build raw SpatialExperiment
################################################################################

sample_info = read_csv(info_path, show_col_types = FALSE)

#   Write transformed image, scalefactors, and coords to a temporary directory
message(Sys.time(), " - Processing ImageJ outputs")
prep_imagej_coords(sample_info, coords_dir)
prep_imagej_image(sample_info, coords_dir)

message(Sys.time(), " - Building raw SPE")
spe = build_spe(sample_info, coords_dir = coords_dir)
colnames(spe) = spe$key

#   Add overlap info, using highest UMI count by capture area to determine
#   which to show (for plots) in case of overlaps
message(Sys.time(), " - Adding overlap info")
spe$sum_umi = Matrix::colSums(assays(spe)$counts)
spe = add_overlap_info(spe, "sum_umi")

################################################################################
#   Filter, log normalize, and save
################################################################################

#   Filter SPE: take only spots in tissue, drop spots with 0 counts for all
#   genes, and drop genes with 0 counts in every spot
spe <- spe[
    rowSums(assays(spe)$counts) > 0,
    (colSums(assays(spe)$counts) > 0) & spe$in_tissue
]

# message(Sys.time(), " - Running quickCluster()")
# spe$scran_quick_cluster <- quickCluster(
#     spe,
#     BPPARAM = MulticoreParam(num_cores),
#     block = spe$capture_area,
#     block.BPPARAM = MulticoreParam(num_cores)
# )
#
# message(Sys.time(), " - Running computeSumFactors()")
# spe <- computeSumFactors(
#     spe,
#    clusters = spe$scran_quick_cluster,
#    BPPARAM = MulticoreParam(num_cores)
#)
##
# message("Quick cluster table:")
# table(spe$scran_quick_cluster)

message("sizeFactors() summary:")
summary(sizeFactors(spe))

message(Sys.time(), " - Running logNormCounts()")
spe <- logNormCounts(spe)

#   Save a copy of the SPE with HDF5-backed assays, which will be important to
#   control memory consumption later#
#message(Sys.time(), " - Saving HDF5-backed object to control memory later")
#spe = saveHDF5SummarizedExperiment(
#    spe, dir = spe_dir, replace = TRUE, as.sparse = TRUE
#)

################################################################################
#   Plot white matter genes to verify smooth transition across capture areas
################################################################################

wm_genes = rownames(spe)[match(wm_genes, rowData(spe)$symbol)]
p = spot_plot(
    spe, sample_id = unique(spe$sample_id)[1], var_name = wm_genes,
    is_discrete = FALSE, assayname = 'logcounts', minCount = 0,
    multi_gene_method = 'pca', spatial = TRUE
)

pdf(file.path(plot_dir, paste0(group,'_WM.pdf')), width = 8, height = 12)
print(p)
dev.off()

session_info()
library(here)
library(ggplot2)
library(SpatialExperiment)
library(SEraster)

# set directories
plots_dir <- here("plots", "Xenium", "07_SEraster")
processed_dir <- here("processed-data", "Xenium", "07_SEraster")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# load xenium data
spe <- readRDS(here("processed-data", "Xenium", "04_dim_reduction", "sce_combined_harmonized_singlecell.rds"))
spe

# specify resolution (set this to whatever you used)
res <- 259

# get unique sample IDs
samples <- unique(spe$brnum)

for (s in samples) {
  message("Processing sample: ", s)
  
  # subset
  spe_sub <- spe[, spe$brnum == s]
  
  # rasterize
  rastGexp <- SEraster::rasterizeGeneExpression(
    spe_sub, 
    assay_name = "counts", 
    resolution = res, 
    square = FALSE, 
    fun = "sum",
    n_threads = 30
  )
  
  # plot
  pdf(here(plots_dir, paste0("SEraster_total_gexp_", s, ".pdf")), width=8, height=6)
  SEraster::plotRaster(rastGexp, name = "total_gexp")
  dev.off()
  
  # save each raster separately
  saveRDS(rastGexp, here(processed_dir, paste0("spe_rasterized_gene_expression_sum_", s, "_", res, "um.rds")))
}
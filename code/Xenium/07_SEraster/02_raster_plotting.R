library(here)
library(ggplot2)
library(SpatialExperiment)
library(SEraster)
library(ggspavis)

processed_dir <- here("processed-data", "Xenium", "07_SEraster")
plots_dir <- here("plots", "Xenium", "07_SEraster")

# load 9017, 9192, 9206, 9280
rast_Br9017 <- readRDS(here(processed_dir, paste0("spe_rasterized_gene_expression_sum_Br9017_100um.rds")))
rast_Br9192 <- readRDS(here(processed_dir, paste0("spe_rasterized_gene_expression_Br9192_100um.rds")))
rast_Br9206 <- readRDS(here(processed_dir, paste0("spe_rasterized_gene_expression_Br9206_100um.rds")))
rast_Br9280 <- readRDS(here(processed_dir, paste0("spe_rasterized_gene_expression_Br9280_100um.rds")))

rast_Br9017$sample_id <- "Br9017"
rast_Br9192$sample_id <- "Br9192"
rast_Br9206$sample_id <- "Br9206"
rast_Br9280$sample_id <- "Br9280"


# cbind
rast_combined <- cbind(rast_Br9017, rast_Br9192, rast_Br9206, rast_Br9280)
rast_combined

pdf(here(plots_dir, "9017_TSHZ1.pdf"), width=10, height=8)
SEraster::plotRaster(rast_Br9017, feature_name = "TSHZ1", name = "TSHZ1")
dev.off()

pdf(here(plots_dir, "9192_TSHZ1.pdf"), width=10, height=8)
SEraster::plotRaster(rast_Br9192, feature_name = "TSHZ1", name = "TSHZ1")
dev.off()

pdf(here(plots_dir, "9206_TSHZ1.pdf"), width=10, height=8)
SEraster::plotRaster(rast_Br9206, feature_name = "TSHZ1", name = "TSHZ1")
dev.off()

pdf(here(plots_dir, "9280_TSHZ1.pdf"), width=10, height=8)   
SEraster::plotRaster(rast_Br9280, feature_name = "TSHZ1", name = "TSHZ1")
dev.off()


# rasters for each samples for CYP26B1, PENK, PDYN, LAMP5, COL25A1
genes_of_interest <- c("CYP26B1", "PENK", "PDYN", "LAMP5", "COL25A1", "MSC", "KLK7", "DLK1", "CD52")

for (gene in genes_of_interest) {
  pdf(here(plots_dir, paste0("9017_", gene, ".pdf")), width=10, height=8)
  p <- SEraster::plotRaster(rast_Br9017, feature_name = gene, name = gene)
  print(p)
  dev.off()
  
  pdf(here(plots_dir, paste0("9192_", gene, ".pdf")), width=10, height=8)
  p <- SEraster::plotRaster(rast_Br9192, feature_name = gene, name = gene)
  print(p)
  dev.off()
  
  pdf(here(plots_dir, paste0("9206_", gene, ".pdf")), width=10, height=8)
  p <- SEraster::plotRaster(rast_Br9206, feature_name = gene, name = gene)
  print(p)
  dev.off()
  
  pdf(here(plots_dir, paste0("9280_", gene, ".pdf")), width=10, height=8)
  p <- SEraster::plotRaster(rast_Br9280, feature_name = gene, name = gene)
  print(p)
  dev.off()
}




# redobut with magma color scale
for (gene in genes_of_interest) {
  pdf(here(plots_dir, paste0("9017_", gene, "_magma.pdf")), width=10, height=8)
  p <- SEraster::plotRaster(rast_Br9017, feature_name = gene, name = gene) +
    scale_fill_viridis_c(option = "magma")
  print(p)
  dev.off() 

    pdf(here(plots_dir, paste0("9192_", gene, "_magma.pdf")), width=10, height=8)
  p <- SEraster::plotRaster(rast_Br9192, feature_name = gene, name = gene) +
    scale_fill_viridis_c(option = "magma")
  print(p)
  dev.off() 

    pdf(here(plots_dir, paste0("9206_", gene, "_magma.pdf")), width=10, height=8)   
    p <- SEraster::plotRaster(rast_Br9206, feature_name = gene, name = gene) +
        scale_fill_viridis_c(option = "magma")
    print(p)
    dev.off()

    pdf(here(plots_dir, paste0("9280_", gene, "_magma.pdf")), width=10, height=8)   
    p <- SEraster::plotRaster(rast_Br9280, feature_name = gene, name = gene) +
        scale_fill_viridis_c(option = "magma")
    print(p)    
    dev.off()
}


gene <- "num_cell"

pdf(here(plots_dir, "9280_total_counts.pdf"), width=10, height=8)
p <- SEraster::plotRaster(rast_Br9280, name = "gene expression")
print(p)
dev.off() 

pdf(here(plots_dir, paste0("9192_", gene, "_magma.pdf")), width=10, height=8)
p <- SEraster::plotRaster(rast_Br9192, feature_name = gene, name = gene) +
    scale_fill_viridis_c(option = "magma")
print(p)
dev.off() 

pdf(here(plots_dir, paste0("9206_", gene, "_magma.pdf")), width=10, height=8)   
p <- SEraster::plotRaster(rast_Br9206, feature_name = gene, name = gene) +
    scale_fill_viridis_c(option = "magma")
print(p)
dev.off()

pdf(here(plots_dir, paste0("9280_", gene, "_magma.pdf")), width=10, height=8)   
p <- SEraster::plotRaster(rast_Br9280, feature_name = gene, name = gene) +
    scale_fill_viridis_c(option = "magma")
print(p)    
dev.off()

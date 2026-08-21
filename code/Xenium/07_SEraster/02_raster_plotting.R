library(here)
library(ggplot2)
library(SpatialExperiment)
library(SEraster)
library(ggspavis)

processed_dir <- here("processed-data", "Xenium", "07_SEraster")
plots_dir <- here("plots", "Xenium", "07_SEraster", "scattermore")

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
# class: SpatialExperiment 
# dim: 366 95158 
# metadata(0):
# assays(1): pixelval
# rownames(366): ABCC9 ADAMTS12 ... ZIC2 ZNF536
# rowData names(0):
# colnames(95158): pixel236 pixel239 ... pixel20778 pixel20779
# colData names(6): num_cell cellID_list ... geometry sample_id
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x y
# imgData names(1): sample_id

pdf(here(plots_dir, "9017_TSHZ1.pdf"), width=10, height=8)
SEraster::plotRaster(rast_Br9017, feature_name = "TSHZ1", name = "TSHZ1") +
  scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()

pdf(here(plots_dir, "9192_TSHZ1.pdf"), width=10, height=8)
SEraster::plotRaster(rast_Br9192, feature_name = "TSHZ1", name = "TSHZ1") +
  scale_fill_gradient(low = "grey", high = "red")
dev.off()

pdf(here(plots_dir, "9206_TSHZ1.pdf"), width=10, height=8)
SEraster::plotRaster(rast_Br9206, feature_name = "TSHZ1", name = "TSHZ1") +
  scale_fill_gradient(low = "white", high = "red") +
  geom_point(stroke = .25, colour="grey")
dev.off()

pdf(here(plots_dir, "9280_SLC17A6.pdf"), width=20, height=18)   
SEraster::plotRaster(rast_Br9280, feature_name = "SLC17A6", name = "SLC17A6") +
  scale_fill_gradient(low = "white", high = "red") +
  geom_point(stroke = .25, colour="grey")

dev.off()


# rasters for each samples for CYP26B1, PENK, PDYN, LAMP5, COL25A1
genes_of_interest <- c("CYP26B1", "PENK", "PDYN", "LAMP5", "COL25A1", "TSHZ1", "PEX5L")

for (gene in genes_of_interest) {
  
  pdf(here(plots_dir, paste0("9192_", gene, ".pdf")), width=20, height=18)
  p <- SEraster::plotRaster(rast_Br9192, feature_name = gene, name = gene)+
  scale_fill_viridis_c(option = "rocket", direction = -1) +
  geom_point(stroke = .25, colour="grey")
  print(p)
  dev.off()
  
  pdf(here(plots_dir, paste0("9280_", gene, ".pdf")), width=20, height=18)
  p <- SEraster::plotRaster(rast_Br9280, feature_name = gene, name = gene)+
  scale_fill_viridis_c(option = "rocket", direction = -1) +
  geom_point(stroke = .25, colour="grey")
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





# ===== scattermoer ====

library(scattermore)
# helper: build a scattermore gene expression plot from a rasterized SPE
plot_gene_scatter <- function(rast_spe, gene, title = gene,
                              pointsize = 6, pixels = c(2048*1.5, 2048*1.5)) {
  plot_df <- data.frame(
    x = spatialCoords(rast_spe)[, 1],
    y = spatialCoords(rast_spe)[, 2],
    expr = assay(rast_spe, "pixelval")[gene, ]
  )

  p <- ggplot(plot_df, aes(x = x, y = y, color = expr)) +
    geom_scattermore(pointsize = pointsize, pixels = pixels) +
    scale_color_viridis_c(option = "rocket", direction = -1) +
    coord_fixed() +
    scale_y_reverse() +
    ggtitle(title) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.title = element_text(size = 12)
    ) +
    labs(color = title)

  p
}

# --- list of samples for easy looping ---
samples <- list(
  Br9017 = rast_Br9017,
  Br9192 = rast_Br9192,
  Br9206 = rast_Br9206,
  Br9280 = rast_Br9280
)

# --- genes of interest ---
genes_of_interest <- c("CYP26B1", "PENK", "PDYN", "LAMP5", "COL25A1", "TSHZ1", "PEX5L")

# reverse rocket plots for all samples x genes
for (sname in names(samples)) {
  for (gene in genes_of_interest) {
    pdf(here(plots_dir, paste0(sname, "_", gene, "_rocket.pdf")), width = 10, height = 8)
    p <- plot_gene_scatter(samples[[sname]], gene, title = gene)
    print(p)
    dev.off()
  }
}

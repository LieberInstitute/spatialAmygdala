suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("RColorBrewer")
    library("ggplot2")
    library("scran")
    library("scater")
    library(patchwork)
})


output_dir <- here("plots", "Xenium", "05_spatial_clustering")

# load
spe <- readRDS(file.path("processed-data", "Xenium", "04_clustering", "Banksy", "Banksy_integrated_res2.0_collapsed_v6.rds"))
spe


# Seperating out Xenium domains

# 1 - WM
# 2* - Glia
# 3 - Glia
# 4 - SST_CHODL
# 5 - Astrocyte
# 6 - Ctx
# 7 - EC?
# 8 - LA
# 9 - HPC?
# 10 - PL
# 11 - BM / CoA
# 12 - BLD
# 13 - small
# 14 - small
# 15 - small
# 16 - MeA/ITC?
# 17 - CeA
# 18 - Endo?
# 19 - Ventricle?
# 20 - small
# 21*
# 22*
# 23* - Astrocyte
# 24 - WEIRD
# 25*

# just get amygdala domains. Ignore weird / uncertain/ endo/ astrocyte/ etc etc

spe$Banksy_domains <- NA

spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 %in% c("Cluster 16", "Cluster 17")] <- "MeA_AI"
spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 == "Cluster 1"] <- "WM"
spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 == "Cluster 6"] <- "Ctx"
spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 == "Cluster 7"] <- "BL"
spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 == "Cluster 8"] <- "LA"
spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 == "Cluster 9"] <- "Sub"
spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 == "Cluster 10"] <- "PL"
spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 == "Cluster 11"] <- "BM_CoA"
spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 == "Cluster 12"] <- "EC"
spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 == "Cluster 17"] <- "CeA"
spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 == "Cluster 18"] <- "Endo"
spe$Banksy_domains[spe$Banksy_res2.0_collapsed_v6 == "Cluster 19"] <- "Ventricle"

# drop NA
spe.na <- spe[, !is.na(spe$Banksy_domains)]

# update new colors for domains
pal_domains <- c(
  "MeA_AI" = "#baf739ff", "BM_CoA" = "#E67E22", "BL" = "#035185ff", 
  "Ctx" = "#d177f5", "EC" = "#7b5289", "LA" = "#F4B400", 
  "Sub" = "#DAB3F5", "PL" = "#f1e438ff", "CeA" = "#197d43ff", 
  "Endo" = "#a3a3a3", "Ventricle" = "#3f3f3f", "WM" = "#BBBBBB")





library(scattermore)

plot_df <- data.frame(
  x = spatialCoords(spe.na)[, 1],
  y = spatialCoords(spe.na)[, 2],
  domain = colData(spe.na)$Banksy_domains,
  sample = colData(spe.na)$brnum
)

pdf_path <- file.path(output_dir, paste0("Banksy_domains_v1.0_scattermore.pdf"))
pdf(pdf_path, width = 10, height = 10)

for (s in unique(plot_df$sample)) {
  p <- ggplot(plot_df[plot_df$sample == s, ], aes(x = x, y = y, color = domain)) +
    geom_scattermore(pointsize = 5, pixels = c(2048*1.5, 2048*1.5)) +
    scale_color_manual(values = pal_domains) +
    coord_fixed() +
    scale_y_reverse() +
    ggtitle(s) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 20),
    ) +
    guides(color = guide_legend(nrow = 2,
                    byrow = TRUE,
                    override.aes = list(size = 5)))
  print(p)
}

dev.off()


# save
saveRDS(spe, file.path("processed-data", "Xenium", "05_spatial_clustering", "Banksy_domains_v1.0.rds"))
library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("bluster")



spe <- readRDS(here("processed-data", "Xenium", "04_clustering","Banksy", "Banksy_integrated_lambda_0.8_res2_gex.rds"))

colnames(colData(spe)) 

# ======= Plotting =======
library(patchwork)


# subset spe objects by brnum
spe.9017 <- spe[,spe$brnum == "Br9017"]
spe.9192 <- spe[,spe$brnum == "Br9192"]
spe.9206 <- spe[,spe$brnum == "Br9206"]
spe.9280 <- spe[,spe$brnum == "Br9280"]


num_groups <- length(unique(spe$clust_HARMONY_M0_lam0.8_k50_res2))

colors <- scCustomize::scCustomize_Palette(
  num_groups,
  ggplot_default_colors = FALSE,
  color_seed = 123
)


pdf(file = here("plots", "Xenium", "05_spatial_clustering", "Banksy_integrated_lambda_0.8_res2.0.pdf"),
    width = 40, height = 10)

p1 <- ggspavis::plotCoords(spe, annotate="clust_HARMONY_M0_lam0.8_k50_res2", in_tissue=NULL, sample_id="brnum") +
    scale_color_manual(values = colors) +
    ggtitle("Banksy Integrated Clustering (lambda=0.8, res=2)") +
    theme(legend.position = "bottom",
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 20),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))

p1
dev.off()




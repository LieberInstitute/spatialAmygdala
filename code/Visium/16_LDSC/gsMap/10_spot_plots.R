#!/usr/bin/env Rscript

library(zellkonverter)
library(SingleCellExperiment)
library(data.table)
library(ggplot2)
library(scales)

proj <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
base_dir <- file.path(proj, "processed-data/Visium/16_LDSC")
h5ad <- file.path(base_dir, "gsMap/ST/Br8325.h5ad")
gsmap_dir <- file.path(base_dir, "gsMap_cond_functional/Br8325/spatial_ldsc")
out_dir <- file.path(base_dir, "gsMap_cond_functional/Br8325/report/spatial_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

traits <- c("SCZ", "MDD", "PTSD", "Height")

spe <- readH5AD(h5ad)
xy <- reducedDim(spe, "spatial")
meta <- data.frame(spot = colnames(spe), x = xy[,1], y = xy[,2])

for(trait in traits){
    f <- file.path(gsmap_dir, paste0("Br8325_", trait, ".csv.gz"))
    d <- fread(f)
    d$neglog10p <- -log10(d$p)
    d$neglog10p[!is.finite(d$neglog10p)] <- 15
    d <- merge(meta, d, by = "spot")
    message(trait, ": ", nrow(d), " spots matched")

    thr <- -log10(0.05)
    p_cap <- 15

    p <- ggplot(d, aes(x = x, y = -y, color = neglog10p)) +
        geom_point(size = 0.3) +
        scale_color_gradientn(
            colors = c("#D2D2D2", "#E4E4E4", "#FFFFFF", "#FDBE85", "#F16913", "#B30000"),
            values = scales::rescale(c(0, thr / 2, thr,
                                    thr + (p_cap - thr) / 3,
                                    thr + 2 * (p_cap - thr) / 3, p_cap)),
            limits = c(0, p_cap), breaks = c(0, thr, 5, 10, 15),
            labels = c("0", "1.3", "5", "10", "≥15"),
            oob = scales::squish, name = expression(-log[10](P))
        ) +
        coord_fixed() +
        theme_void() +
        ggtitle(trait) +
        theme(plot.title = element_text(hjust = 0.5, size = 14))

    ggsave(file.path(out_dir, paste0("Br8325_", trait, "_spatial.pdf")), p, width = 5, height = 5)
}
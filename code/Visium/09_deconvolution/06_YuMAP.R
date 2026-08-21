library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialFeatureExperiment)
library(scater)
library(scran)
library(Voyager)
library(ggspavis)
library(harmony)
library(escheR)
library(spacexr)

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")

# load
load(here("processed-data", "snRNAseq", "Yu_et_al", "yu_sce_gtf.rda"))
sce.amy

counts(sce.amy) <- round(expm1(counts(sce.amy)))

sce <- sce.amy

library(SingleCellExperiment)
library(ggplot2)
library(scattermore)

# ============================================================
# Assign each inhibitory ident to the annotation-bar subclasses
# ============================================================

cge_idents <- c(
    "Human_CALCR LHX8",
    "Human_LAMP5 COL14A1",
    "Human_LAMP5 NDNF",
    "Human_VIP ABI3BP",
    "Human_VIP NDNF"
)

mge_idents <- c(
    "Human_PVALB ADAMTS5",
    "Human_SST EPYC",
    "Human_SST HGF"
)

lge_idents <- c(
    "Human_DRD2 ISL1",
    "Human_DRD2 PAX6",
    "Human_HTR3A DRD2",
    "Human_PRKCD",
    "Human_TSHZ1 CALCRL",
    "Human_TSHZ1 SEMA3C"
)

ident_chr <- as.character(sce$ident)
celltype_chr <- as.character(sce$celltype)
all_idents <- unique(ident_chr)

# Determine the broad plotting group for each ident
ident_group <- setNames(rep(NA_character_, length(all_idents)), all_idents)

ident_group[all_idents %in% unique(ident_chr[celltype_chr == "ExN"])] <-
    "ExN"

ident_group[all_idents %in% cge_idents] <-
    "CGE InN"

ident_group[all_idents %in% mge_idents] <-
    "MGE InN"

ident_group[all_idents %in% lge_idents] <-
    "LGE InN"

ident_group[all_idents %in% unique(
    ident_chr[!celltype_chr %in% c("ExN", "InN")]
)] <- "Non-neuron"

# Stop if any ident was not assigned
unassigned <- names(ident_group)[is.na(ident_group)]

if (length(unassigned) > 0) {
    stop(
        "The following idents were not assigned:\n",
        paste(unassigned, collapse = "\n")
    )
}

# ============================================================
# Generate varying shades within each annotation-bar color
# ============================================================

make_shades <- function(group, light, dark) {
    ids <- names(ident_group)[ident_group == group]

    setNames(
        grDevices::colorRampPalette(c(light, dark))(length(ids)),
        ids
    )
}

ident_cols <- c(
    make_shades(
        "ExN",
        light = "#FDBE85",
        dark  = "#D95F0E"
    ),
    make_shades(
        "CGE InN",
        light = "#FCAE91",
        dark  = "#CB181D"
    ),
    make_shades(
        "MGE InN",
        light = "#BDD7E7",
        dark  = "#08519C"
    ),
    make_shades(
        "LGE InN",
        light = "#C7E9C0",
        dark  = "#238B45"
    ),
    make_shades(
        "Non-neuron",
        light = "#E6E1F2",
        dark  = "#9E8AC7"
    )
)

# ============================================================
# Extract UMAP coordinates and plot with scattermore
# ============================================================

umap <- reducedDim(sce, "UMAP")

umap_df <- data.frame(
    UMAP1 = umap[, 1],
    UMAP2 = umap[, 2],
    ident = ident_chr
)

p <- ggplot(
    umap_df,
    aes(
        x = UMAP1,
        y = UMAP2,
        colour = ident
    )
) +
    scattermore::geom_scattermore(
        pointsize = 4,
        pixels = c(2000, 2000)
    ) +
    scale_colour_manual(
        values = ident_cols
    ) +
    coord_equal() +
    theme_classic() +
    theme(
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

ggsave(
    filename = file.path(plots_dir, "UMAP_Yu_idents.pdf"),
    plot = p,
    width = 7,
    height = 7,
    device = cairo_pdf
)
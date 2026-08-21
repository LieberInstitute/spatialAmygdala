library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialFeatureExperiment)
library(spacexr)

processed_dir <- here("processed-data", "Xenium", "06_label_transfer")
plot_dir <- here("plots", "Xenium", "06_label_transfer")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# Load reference to build ident -> celltype mapping
# --------------------------------------------------------------------------
load(here("processed-data", "snRNAseq", "Yu_et_al", "yu_sce_gtf.rda"))
sce <- sce.amy

## Build a lookup table: ident -> broad celltype
ident_to_celltype <- data.frame(
    ident = sce$ident,
    celltype = sce$celltype,
    stringsAsFactors = FALSE
) %>%
    distinct()

rm(sce.amy, sce)

# --------------------------------------------------------------------------
# Load xenium data
# --------------------------------------------------------------------------
spe <- readRDS(here("processed-data", "Xenium", "05_spatial_clustering",
                     "Banksy_domains_v1.0.rds"))



# --------------------------------------------------------------------------
# Load RCTD results and combine
# --------------------------------------------------------------------------
sample_ids <- c("Br9017", "Br9192", "Br9206", "Br9280")

rctd_list <- lapply(sample_ids, function(s) {
    readRDS(here(processed_dir, paste0("rctd_results_", s, ".rds")))
})
names(rctd_list) <- sample_ids


# --------------------------------------------------------------------------
# Add RCTD labels to spe — match by sample + original barcode
# --------------------------------------------------------------------------
sample_ids <- sort(unique(spe$brnum))

## For each sample, match RCTD barcodes to spe cells of that sample
spe$rctd_ident <- NA_character_
spe$rctd_celltype <- NA_character_

for (s in sample_ids) {
    ## Cells in spe belonging to this sample
    spe_idx <- which(spe$brnum == s)
    
    ## Get RCTD results for this sample
    res_df <- rctd_list[[s]]@results$results_df
    rctd_barcodes <- rownames(res_df)
    
    ## The original (pre-make.unique) barcodes for this sample's cells
    ## Since subsetting by sample should give unique barcodes within sample,
    ## strip any make.unique suffixes or just re-derive from spatialCoords
    original_barcodes <- colnames(spe)[spe_idx]
    
    ## Try direct match first
    m <- match(original_barcodes, rctd_barcodes)
    
    matched <- !is.na(m)
    message(s, ": matched ", sum(matched), " / ", length(spe_idx), " cells")
    
    spe$rctd_ident[spe_idx[matched]] <- as.character(res_df$first_type[m[matched]])
}

## Now add broad celltype
spe$rctd_celltype <- ident_to_celltype$celltype[
    match(spe$rctd_ident, ident_to_celltype$ident)
]


# save
saveRDS(spe, here(processed_dir, "spe_with_rctd_labels_ITCmn.rds"))


# --------------------------------------------------------------------------
# Identify excitatory / inhibitory idents
# --------------------------------------------------------------------------
# Drop "Human_" prefix from idents for cleaner plotting
spe$rctd_ident <- gsub("^Human_", "", spe$rctd_ident)
ident_to_celltype$ident <- gsub("^Human_", "", ident_to_celltype$ident)

exn_idents <- ident_to_celltype %>% filter(celltype == "ExN") %>% pull(ident)
inn_idents <- ident_to_celltype %>% filter(celltype == "InN") %>% pull(ident)

# add ITC_1, ITC_2, ITC_3 to inhibitory list
#inn_idents <- c(inn_idents, "ITC_1", "ITC_2", "ITC_3")

message("Excitatory subtypes: ", paste(exn_idents, collapse = ", "))
message("Inhibitory subtypes: ", paste(inn_idents, collapse = ", "))



# --------------------------------------------------------------------------
# Plotting function
# --------------------------------------------------------------------------
plot_rctd_spatial <- function(spe, color_col, title_suffix = "",
                               subset_idents = NULL, palette = NULL) {
    plots <- list()
    for (s in sample_ids) {
        spe_sub <- spe[, spe$brnum == s]

        ## Subset to specific idents if requested
        if (!is.null(subset_idents)) {
            spe_sub <- spe_sub[, !is.na(spe_sub$rctd_ident) &
                                   spe_sub$rctd_ident %in% subset_idents]
        }

        ## Remove cells with no RCTD assignment
        spe_sub <- spe_sub[, !is.na(colData(spe_sub)[[color_col]])]

        coords <- as.data.frame(spatialCoords(spe_sub))
        colnames(coords) <- c("x", "y")
        coords[[color_col]] <- colData(spe_sub)[[color_col]]

        p <- ggplot(coords, aes(x = x, y = y, color = .data[[color_col]])) +
            geom_point(size = 1, stroke = 0) +
            coord_equal() +
            theme_void(base_size = 14) +
            theme(
                legend.position = "right",
                legend.text = element_text(size = 8),
                plot.title = element_text(face = "bold")
            ) +
            guides(color = guide_legend(override.aes = list(size = 3))) +
            ggtitle(paste0(s, " ", title_suffix))

        ## Apply a custom color palette if requested
        if (!is.null(palette)) {
            p <- p + scale_color_manual(values = palette)
        }

        plots[[s]] <- p
    }
    plots
}

# custom color palette for ITC subtypes
pal <- c("ITC_1" = "#fd0d00",
         "ITC_2" = "#f5b6b3",
         "ITC_3" = "#5d0500")

# --------------------------------------------------------------------------
# PDF 1: All cell types (rctd_ident)
# --------------------------------------------------------------------------
message("Plotting all cell types...")
pdf(here(plot_dir, "RCTD_Yu_all_celltypes.pdf"), width = 14, height = 10)
plots_all <- plot_rctd_spatial(spe, "rctd_ident", "- All Cell Types")
for (p in plots_all) print(p)
dev.off()

# --------------------------------------------------------------------------
# PDF 2: Excitatory neurons only
# --------------------------------------------------------------------------
message("Plotting excitatory neurons...")
pdf(here(plot_dir, "RCTD_Yu_excitatory_neurons.pdf"), width = 14, height = 10)
plots_exn <- plot_rctd_spatial(spe, "rctd_ident", "- Excitatory Neurons",
                                subset_idents = exn_idents)
for (p in plots_exn) print(p)
dev.off()

# --------------------------------------------------------------------------
# PDF 3: Inhibitory neurons only
# --------------------------------------------------------------------------
message("Plotting inhibitory neurons...")
pdf(here(plot_dir, "RCTD_Yu_inhibitory_neurons_ITCmn.pdf"), width = 10, height = 7)
plots_inn <- plot_rctd_spatial(spe, "rctd_ident", "- Inhibitory Neurons",
                                subset_idents = inn_idents)
for (p in plots_inn) print(p)
dev.off()

# --------------------------------------------------------------------------
# PDF 3b: ITC subtypes only (custom palette)
# --------------------------------------------------------------------------
message("Plotting ITC subtypes...")
pdf(here(plot_dir, "RCTD_Yu_ITC_subtypes.pdf"), width = 10, height = 7)
plots_itc <- plot_rctd_spatial(spe, "rctd_ident", "- ITC Subtypes",
                                subset_idents = c("ITC_1", "ITC_2", "ITC_3"),
                                palette = pal)
for (p in plots_itc) print(p)
dev.off()

message("Done! PDFs saved to: ", plot_dir)


# --------------------------------------------------------------------------
# PDF 4: All neurons (excitatory + inhibitory)
# --------------------------------------------------------------------------
message("Plotting all neurons...")
all_neuron_idents <- c(exn_idents, inn_idents)
pdf(here(plot_dir, "RCTD_all_neurons.pdf"), width = 14*0.75, height = 10*0.75)
plots_neurons <- plot_rctd_spatial(spe, "rctd_ident", "- All Neurons",
                                    subset_idents = all_neuron_idents)
for (p in plots_neurons) print(p)
dev.off()


# --------------------------------------------------------------------------
# Br9280 all neurons — rotated 180° + mirrored vertically
# --------------------------------------------------------------------------
spe_sub <- spe[, spe$brnum == "Br9280"]
neuron_idents <- c(exn_idents, inn_idents)
spe_sub <- spe_sub[, !is.na(spe_sub$rctd_ident) &
                       spe_sub$rctd_ident %in% neuron_idents]

coords <- as.data.frame(spatialCoords(spe_sub))
colnames(coords) <- c("x", "y")
coords$rctd_ident <- spe_sub$rctd_ident

## Rotate 180°: negate x and y
## Mirror vertically: negate y again
## Net effect: negate x only
coords$y <- -coords$y

# --------------------------------------------------------------------------
# PDF 5: Selected neuron subtypes (custom palette)
# --------------------------------------------------------------------------
selected_pal <- c(
    "VGLL3 MEPE"    = "#E8780C",  # orange
    "LAMP5 COL25A1" = "#E6E64D",  # yellow
    "LAMP5 ABO"     = "#1F5C99",  # dark blue
    "HGF NPSR1"     = "#E89A3C",  # light orange-tan
    "HGF ESR1"      = "#F5C77E",  # pale orange
    "PRKCD"         = "#1A7A3C",  # dark green
    "STRIP2"        = "#7AC043",  # medium green
    "SATB2 CALCRL"   = "#54a6fd",  # light green/lime
    "TFAP2C"        = "#C5E063",  # light green/lime
    "TSHZ1 SEMA3C"  = "#7A0A0A",  # dark red/maroon
    "TSHZ1 CALCRL"   = "#F07060",  # salmon/coral
    "Oligo_3 OPALIN" = "#B0B0B0"   # grey
)

selected_idents <- names(selected_pal)

message("Plotting selected neuron subtypes...")
pdf(here(plot_dir, "RCTD_selected_neurons.pdf"), width = 14*0.75, height = 10*0.75)
plots_selected <- plot_rctd_spatial(spe, "rctd_ident", "- Selected Neurons",
                                     subset_idents = selected_idents,
                                     palette = selected_pal)
for (p in plots_selected) print(p)
dev.off()






# ======== HEATMAP ========


suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
})

# --------------------------------------------------------------------------
# Strip "Human_" prefix from RCTD ident labels
# --------------------------------------------------------------------------
spe$rctd_ident <- gsub("^Human_", "", spe$rctd_ident)
ident_to_celltype$ident <- gsub("^Human_", "", ident_to_celltype$ident)

spe$rctd_celltype <- ident_to_celltype$celltype[
    match(spe$rctd_ident, ident_to_celltype$ident)
]

# --------------------------------------------------------------------------
# Build RCTD_broad_celltype: split InN into CGE/MGE/LGE based on fine labels
# --------------------------------------------------------------------------
ident_to_celltype$RCTD_broad <- as.character(ident_to_celltype$celltype)

## MGE: SST, PVALB
ident_to_celltype$RCTD_broad[
    grepl("SST|PVALB", ident_to_celltype$ident)
] <- "MGE"

## CGE: LAMP5, VIP, CCK, CALCR
ident_to_celltype$RCTD_broad[
    grepl("LAMP5|VIP|CCK|CALCR", ident_to_celltype$ident)
] <- "CGE"

## LGE: remaining inhibitory neurons not in MGE/CGE
ident_to_celltype$RCTD_broad[
    ident_to_celltype$celltype == "InN" &
        !grepl("SST|PVALB|LAMP5|VIP|CCK|LHX8", ident_to_celltype$ident)
] <- "LGE"

## Excitatory and Non-neuronal
ident_to_celltype$RCTD_broad[ident_to_celltype$celltype == "ExN"] <- "Excitatory"
ident_to_celltype$RCTD_broad[
    !ident_to_celltype$RCTD_broad %in% c("Excitatory", "CGE", "MGE", "LGE")
] <- "Non-neuronal"

# --------------------------------------------------------------------------
# Compute enrichment matrix: domain x fine cell type
# --------------------------------------------------------------------------
has_rctd <- !is.na(spe$rctd_ident)
df <- data.frame(
    domain = as.character(spe$Banksy_domains[has_rctd]),
    celltype = spe$rctd_ident[has_rctd],
    stringsAsFactors = FALSE
)

obs <- table(df$domain, df$celltype)
expected <- outer(rowSums(obs), colSums(obs)) / sum(obs)
log2_enrich <- log2((obs + 1) / (expected + 1))

## Z-score per cell type (scale columns)
enrich_mat <- base::scale(as.matrix(log2_enrich))

# --------------------------------------------------------------------------
# Map fine cell types -> broad classes and reorder columns
# --------------------------------------------------------------------------
rctd_celltypes <- colnames(enrich_mat)
broad_labels <- ident_to_celltype$RCTD_broad[
    match(rctd_celltypes, ident_to_celltype$ident)
]

keep_cols <- !is.na(broad_labels)
enrich_mat <- enrich_mat[, keep_cols, drop = FALSE]
rctd_celltypes <- colnames(enrich_mat)
broad_labels <- broad_labels[keep_cols]

annotation_col <- data.frame(BroadClass = broad_labels, stringsAsFactors = FALSE)
rownames(annotation_col) <- rctd_celltypes

## Reorder columns by broad class
desired_order <- c("Excitatory", "CGE", "MGE", "LGE", "Non-neuronal")
ordered_celltypes <- unlist(lapply(desired_order, function(cls) {
    sort(rownames(annotation_col)[annotation_col$BroadClass == cls])
}))

enrich_mat_ordered <- enrich_mat[, ordered_celltypes, drop = FALSE]
annotation_col_ordered <- annotation_col[ordered_celltypes, , drop = FALSE]

# --------------------------------------------------------------------------
# Reorder rows (domains)
# -------------------------------------------------------------------------
domain_order <- c("LA", "BL", "PL", "EC", "Sub", "Ctx",
                   "MeA_AI", "BM_CoA", "CeA", "Endo", "Ventricle", "WM")
domain_order_present <- intersect(domain_order, rownames(enrich_mat_ordered))
enrich_mat_ordered <- enrich_mat_ordered[domain_order_present, , drop = FALSE]

annotation_row <- data.frame(
    Domain = factor(domain_order_present, levels = domain_order_present),
    stringsAsFactors = FALSE
)
rownames(annotation_row) <- domain_order_present

# --------------------------------------------------------------------------
# Colors
# --------------------------------------------------------------------------
anno_colors <- list(
    BroadClass = c(
        "CGE" = "#E15759",
        "MGE" = "#4E79A7",
        "LGE" = "#59A14F",
        "Excitatory" = "#F28E2B",
        "Non-neuronal" = "#B07AA1"
    )
)

domain_pal <- c(
    BM_CoA    = "#E67E22",
    PL        = "#f1e438",
    BL        = "#035185",
    LA        = "#F4B400",
    Sub       = "#5DA5DA",
    EC        = "#9B59B6",
    Ctx       = "#d6a8f8",
    MeA_AI    = "#baf739",
    CeA       = "#197d43",
    Endo      = "#444444",
    Ventricle = "#A0522D",
    WM        = "#BBBBBB"
)

# --------------------------------------------------------------------------
# Build annotations + heatmap
# --------------------------------------------------------------------------
col_ha <- HeatmapAnnotation(
    BroadClass = factor(annotation_col_ordered$BroadClass, levels = desired_order),
    which = "column",
    col = list(BroadClass = anno_colors$BroadClass),
    show_annotation_name = FALSE
)

row_ha <- rowAnnotation(
    Domain = annotation_row$Domain,
    col = list(Domain = domain_pal),
    show_annotation_name = FALSE
)

rng <- range(enrich_mat_ordered, na.rm = TRUE, finite = TRUE)
heatmap_colors <- colorRamp2(
    breaks = c(rng[1], 0, rng[2]),
    colors = c("purple", "white", "darkgreen")
)

ht <- Heatmap(
    enrich_mat_ordered,
    name = "Z-score",
    col = heatmap_colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    bottom_annotation = col_ha,
    right_annotation = row_ha,
    column_split = factor(annotation_col_ordered$BroadClass, levels = desired_order),
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(title = "Z-score", legend_direction = "vertical"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        v <- enrich_mat_ordered[i, j]
        if (is.finite(v) && v > 2) {
            grid::grid.text("X", x = x, y = y,
                             gp = grid::gpar(fontsize = 10, fontface = "bold", col = "white"))
        }
    }
)

pdf(here(plot_dir, "RCTD_celltype_enrichment_by_domain.pdf"),
    width = 12, height = 6)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("Enrichment heatmap saved to: ", plot_dir)
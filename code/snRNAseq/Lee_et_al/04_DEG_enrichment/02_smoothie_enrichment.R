library(here)
library(dplyr)
library(tidyr)
library(spatialLIBD)   # for layer_matrix_plot styling

# ==============================================================================
# 1. Smoothie modules (gene -> module)
# ==============================================================================
results_dir <- here("processed-data","Visium","14_smoothie_coexpression","results")
modules_df  <- read.csv(file.path(results_dir, "modules_df_0.8_9.csv"))
# columns: name (gene symbol), module_label (int)

# Module -> domain label map (from your network plotting script)
mod_domain <- c(
  "6"="LA", "27"="BLD", "5"="PL", "10"="CoA",
  "13"="CeA", "19"="MeA", "20"="AI", "16"="HPC"
)

# Build named list of module gene sets — ONLY the annotated domain-specific modules
module_sets <- split(modules_df$name, modules_df$module_label)
names(module_sets) <- paste0("M", names(module_sets))
module_sets <- module_sets[paste0("M", names(mod_domain))]   # keep domain modules only
# Friendly labels
annot_names <- setNames(paste0("M", names(mod_domain), "_", mod_domain),
                        paste0("M", names(mod_domain)))

# Universe = all genes assigned to any module (Smoothie background)
universe <- unique(modules_df$name)

# ==============================================================================
# 2. Girgenti DEGs (symbol-keyed)
# ==============================================================================
degs <- read.csv(here("processed-data","snRNAseq","Lee_at_al",
                      "Supplementary Data","Supplementary Data 5.csv"))
gene_col <- grep("gene|symbol", colnames(degs), value=TRUE, ignore.case=TRUE)[1]
ct_col   <- grep("cell|cluster|type", colnames(degs), value=TRUE, ignore.case=TRUE)[1]
fc_col   <- grep("log2|logfc|fc|coef", colnames(degs), value=TRUE, ignore.case=TRUE)[1]

degs <- degs |>
  mutate(direction = ifelse(.data[[fc_col]] > 0, "up", "down"),
         set_name  = paste0(.data[[ct_col]], "_", direction))

deg_sets <- lapply(split(degs[[gene_col]], degs$set_name), unique)
deg_sets <- deg_sets[lengths(deg_sets) >= 10]   # drop underpowered sets

# ==============================================================================
# 3. Fisher's exact: each DEG set x each module
# ==============================================================================
fisher_enrich <- function(deg_genes, mod_genes, universe) {
  deg_genes <- intersect(deg_genes, universe)
  mod_genes <- intersect(mod_genes, universe)
  a <- length(intersect(deg_genes, mod_genes))         # in both
  b <- length(setdiff(deg_genes, mod_genes))           # DEG not module
  c <- length(setdiff(mod_genes, deg_genes))           # module not DEG
  d <- length(universe) - a - b - c                    # neither
  ft <- fisher.test(matrix(c(a,b,c,d), 2), alternative = "greater")
  c(OR = unname(ft$estimate), p = ft$p.value, overlap = a)
}

res <- expand.grid(deg = names(deg_sets), mod = names(module_sets),
                   stringsAsFactors = FALSE)
res <- cbind(res, t(mapply(function(dg, md)
  fisher_enrich(deg_sets[[dg]], module_sets[[md]], universe),
  res$deg, res$mod)))
res$fdr <- p.adjust(res$p, method = "BH")

write.csv(res, here("processed-data","snRNAseq","dreamlet",
                    "girgenti_smoothie_module_enrichment.csv"), row.names = FALSE)

# ==============================================================================
# 4. LIBD-style heatmap via ComplexHeatmap (modules x DEG sets), per direction
# ==============================================================================
library(ComplexHeatmap)
library(circlize)
library(grid)

# Domain color palette (for row annotation)
pal <- c(
  AI = "#D62728", BM = "#E67E22", BLD = "#9B59B6", PL = "#f1e438ff",
  BL = "#035185ff", LA = "#F4B400", CoA = "#5DA5DA", CeA = "#197d43ff",
  MeA = "#baf739ff", HPC = "#d6a8f8ff", CHAT = "#A0522D",
  Endothelial = "#444444", WM.1 = "#BBBBBB", WM.2 = "#DDDDDD", CLA = "#FF69B4"
)

make_heatmap <- function(direction, high_col, fname) {
  dsub <- res |> filter(grepl(paste0("_", direction, "$"), deg))
  dsub$deg_lab <- gsub(paste0("_", direction, "$"), "", dsub$deg)
  dsub$mod_lab <- ifelse(dsub$mod %in% names(annot_names),
                         annot_names[dsub$mod], dsub$mod)

  pmat <- dsub |>
    mutate(neglogp = -log10(p)) |>
    select(mod_lab, deg_lab, neglogp) |>
    pivot_wider(names_from = mod_lab, values_from = neglogp) |>
    as.data.frame()
  rownames(pmat) <- pmat$deg_lab; pmat$deg_lab <- NULL
  pmat <- t(as.matrix(pmat))

  ormat <- dsub |>
    mutate(orlab = ifelse(p < 0.05 & OR > 1, sprintf("%.1f", OR), "")) |>
    select(mod_lab, deg_lab, orlab) |>
    pivot_wider(names_from = mod_lab, values_from = orlab) |>
    as.data.frame()
  rownames(ormat) <- ormat$deg_lab; ormat$deg_lab <- NULL
  ormat <- t(as.matrix(ormat))[rownames(pmat), colnames(pmat)]

  set_sizes <- lengths(deg_sets)[paste0(colnames(pmat), "_", direction)]
  top_anno <- HeatmapAnnotation(
    SetSize = anno_barplot(set_sizes, gp = gpar(fill = "grey40"),
                           height = unit(1.5, "cm")),
    annotation_name_side = "left"
  )

  # Right-side row annotation: domain color per module (parse domain from M#_DOMAIN)
  row_domains <- sub("^M\\d+_", "", rownames(pmat))
  row_anno <- rowAnnotation(
    Domain = row_domains,
    col = list(Domain = pal[intersect(names(pal), unique(row_domains))]),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  # Shared color scale across up/down panels — cap at PCAP; values above -> max color
  PCAP <- 6
  pmat_capped <- pmin(pmat, PCAP)
  ramp_pal  <- if (direction == "up") "YlOrRd" else "Blues"
  libd_cols <- c("white",
                 grDevices::colorRampPalette(
                   RColorBrewer::brewer.pal(9, ramp_pal))(50))
  col_fun <- colorRamp2(seq(0, PCAP, length.out = length(libd_cols)), libd_cols)

  ht <- Heatmap(
    pmat_capped,
    name = "-log10(p-val)",
    col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    rect_gp = gpar(col = "black", lwd = 0.6),
    top_annotation = top_anno,
    right_annotation = row_anno,
    column_names_rot = 45,
    row_names_side = "left",
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (ormat[i, j] != "")
        grid.text(ormat[i, j], x, y, gp = gpar(fontsize = 9))
    },
    heatmap_legend_param = list(at = c(0, 2, 4, 6),
                                legend_height = unit(3, "cm"))
  )

  pdf(file.path(here("plots","snRNAseq","girgenti"), fname),
      width = 8, height = 6)
  draw(ht, column_title = paste0("Girgenti AUD DEGs (", direction,
                                 ") vs Smoothie domain modules"))
  dev.off()
}

dir.create(here("plots","snRNAseq","girgenti"), recursive=TRUE, showWarnings=FALSE)
make_heatmap("up",   "#b2182b", "girgenti_smoothie_enrichment_UP.pdf")
make_heatmap("down", "#08519c", "girgenti_smoothie_enrichment_DOWN.pdf")
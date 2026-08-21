library(here)
library(spatialLIBD)
library(SingleCellExperiment)
library(dplyr)

# ===== 1. Spatial modeling results (BS k=16) =====
load(here("processed-data","Visium","08_marker_genes",
          "BS_k16_modeling_results.Rdata"))
print(names(modeling_results))
print(colnames(modeling_results$enrichment))

# ===== 2. symbol -> ensembl map from Girgenti rowData =====
sce <- readRDS(here("processed-data","snRNAseq","Lee_at_al","girgenti_sce.rds"))
sym2ens <- setNames(rowData(sce)$featureid, rownames(sce))

# ===== 3. Published Girgenti DEGs =====
degs <- read.csv(here("processed-data","snRNAseq","Lee_at_al",
                      "Supplementary Data","Supplementary Data 5.csv"))
print(colnames(degs)); print(head(degs))

gene_col <- grep("gene|symbol", colnames(degs), value = TRUE, ignore.case = TRUE)[1]
ct_col   <- grep("cell|cluster|type", colnames(degs), value = TRUE, ignore.case = TRUE)[1]
fc_col   <- grep("log2|logfc|fc|coef", colnames(degs), value = TRUE, ignore.case = TRUE)[1]

# ===== 4. Build gene sets: (cell type x direction) -> ensembl =====
degs <- degs |>
  mutate(direction = ifelse(.data[[fc_col]] > 0, "up", "down"),
         set_name  = paste0(.data[[ct_col]], "_", direction),
         ensembl   = sym2ens[.data[[gene_col]]]) |>
  filter(!is.na(ensembl))

gene_list <- lapply(split(degs$ensembl, degs$set_name), unique)

# ===== 5. Fisher's exact enrichment, split by direction =====
up_sets   <- gene_list[grepl("_up$",   names(gene_list))]
down_sets <- gene_list[grepl("_down$", names(gene_list))]

enrich_up <- gene_set_enrichment(
  gene_list = up_sets, modeling_results = modeling_results,
  model_type = "enrichment", fdr_cut = 0.05)
enrich_down <- gene_set_enrichment(
  gene_list = down_sets, modeling_results = modeling_results,
  model_type = "enrichment", fdr_cut = 0.05)

# ===== 6. Plot — UP (default red) and DOWN (blue) =====
out_plots <- here("plots","snRNAseq","girgenti")
dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)

# Palettes: UP = default YlOrRd (red), DOWN = blues
red_pal  <- c("white", grDevices::colorRampPalette(
                RColorBrewer::brewer.pal(9, "YlOrRd"))(50))
blue_pal <- c("white", grDevices::colorRampPalette(
                RColorBrewer::brewer.pal(9, "Blues"))(50))

pdf(file.path(out_plots, "girgenti_AUD_enrichment_BSk16_UP.pdf"),
    width = 10, height = 8)
gene_set_enrichment_plot(enrich_up, PThresh = 12, ORcut = 2,
                         mypal = red_pal,
                         xlabs = gsub("_up$", "", unique(enrich_up$ID)))
dev.off()

pdf(file.path(out_plots, "girgenti_AUD_enrichment_BSk16_DOWN.pdf"),
    width = 10, height = 8)
gene_set_enrichment_plot(enrich_down, PThresh = 12, ORcut = 2,
                         mypal = blue_pal,
                         xlabs = gsub("_down$", "", unique(enrich_down$ID)))
dev.off()

write.csv(rbind(transform(enrich_up, dir = "up"),
                transform(enrich_down, dir = "down")),
          here("processed-data","snRNAseq","dreamlet",
               "girgenti_spatial_enrichment_BSk16.csv"),
          row.names = FALSE)
library(here)
library(SingleCellExperiment)
library(ggplot2)
library(ggsankey)
library(dplyr)
library(RColorBrewer)

sce <- readRDS(here("processed-data","snRNAseq","Lee_at_al","girgenti_INH_yu.rds"))

df <- data.frame(
  girgenti_major   = sce$celltype,
  girgenti_subtype = sce$subtype,
  yu_fine          = gsub("^Human_", "", sce$yu_inh_predicted.id)
)

make_sankey <- function(major_label, title) {
  d <- df |>
    filter(girgenti_major == major_label) |>
    count(girgenti_subtype, yu_fine, name = "value") |>
    group_by(girgenti_subtype) |>
    mutate(prop = value / sum(value)) |>
    ungroup() |>
    filter(prop > 0.02)

  # Palette keyed to Girgenti subtypes
  subs <- sort(unique(d$girgenti_subtype))
  pal  <- setNames(brewer.pal(max(3, length(subs)), "Set2")[seq_along(subs)],
                   subs)

  # Color each Yu node by its DOMINANT Girgenti source (largest inflow)
  yu_dominant <- d |>
    group_by(yu_fine) |>
    slice_max(value, n = 1, with_ties = FALSE) |>
    transmute(yu_fine, src = girgenti_subtype) |>
    ungroup()
  yu_pal <- setNames(pal[yu_dominant$src], yu_dominant$yu_fine)

  full_pal <- c(pal, yu_pal)

  d_long <- bind_rows(
    d |> transmute(x = "Girgenti", node = girgenti_subtype,
                   next_x = "Yu",   next_node = yu_fine, value = prop),
    d |> transmute(x = "Yu",       node = yu_fine,
                   next_x = NA,    next_node = NA,       value = prop)
  )
  d_long$x      <- factor(d_long$x,      levels = c("Girgenti", "Yu"))
  d_long$next_x <- factor(d_long$next_x, levels = c("Girgenti", "Yu"))

  n_left  <- length(unique(d$girgenti_subtype))
  n_right <- length(unique(d$yu_fine))

  ggplot(d_long,
         aes(x = x, next_x = next_x, node = node, next_node = next_node,
             value = value, fill = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = "grey40",
                space = 0.08, width = 0.12) +
    geom_sankey_text(aes(label = node), size = 3.2,
                     hjust = c(rep(1.1, n_left), rep(-0.1, n_right)),
                     space = 0.08) +
    scale_fill_manual(values = full_pal) +
    labs(title = title) +
    theme_sankey(base_size = 12) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
}

p_inh <- make_sankey("INH", "Inhibitory: Girgenti -> Yu")

out_plots <- here("plots","snRNAseq","girgenti")
dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_plots, "sankey_INH.pdf"), p_inh, width = 9, height = 7)
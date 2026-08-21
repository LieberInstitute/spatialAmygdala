library(here)
library(igraph)
library(ggraph)
library(graphlayouts)
library(tidygraph)
library(ggforce)
library(dplyr)

# ==============================================================================
# Configuration
# ==============================================================================
results_dir <- here("processed-data", "Visium", "14_smoothie_coexpression", "results")
plots_dir <- here("plots", "Visium", "14_smoothie_coexpression")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

network_df <- read.csv(file.path(results_dir, "network", "network_PCC0.8_clustPow9.csv"))
modules_df <- read.csv(file.path(results_dir, "modules_df_0.8_9.csv"))


# ==============================================================================
# Domain palette
# ==============================================================================
pal <- c(
    AI          = "#D62728",
    BM          = "#E67E22",
    BLD         = "#9B59B6",
    PL          = "#f1e438ff",
    BL          = "#035185ff",
    LA          = "#F4B400",
    CoA         = "#5DA5DA",
    CeA         = "#197d43ff",
    MeA         = "#baf739ff",
    HPC         = "#d6a8f8ff",
    CHAT        = "#A0522D",
    Endothelial = "#444444",
    WM.1        = "#BBBBBB",
    WM.2        = "#DDDDDD",
    CLA         = "#FF69B4"
)

# ==============================================================================
# Manual module selection + domain assignment
# ==============================================================================
keep_modules <- c(7, 5, 10, 16, 27, 6, 13, 19, 20, 15)

domain_lookup <- c(
    "7"  = "CLA",
    "5"  = "PL",
    "10" = "CoA",
    "16" = "HPC",
    "27" = "BLD",
    "6"  = "LA",
    "13" = "CeA",
    "19" = "MeA",
    "20" = "AI",
    "15" = "CHAT"
)

kept_genes <- modules_df %>%
    filter(module_label %in% keep_modules) %>%
    pull(name)

mod_sizes <- modules_df %>%
    filter(module_label %in% keep_modules) %>%
    group_by(module_label) %>%
    summarise(n_genes = n(), .groups = "drop")

# Hull labels
hull_labels <- mod_sizes %>%
    mutate(label = paste0("M", module_label, " - ",
                          domain_lookup[as.character(module_label)],
                          " (", n_genes, ")"))
label_lookup <- setNames(hull_labels$label, as.character(hull_labels$module_label))

# ==============================================================================
# Build network (with light pruning)
# ==============================================================================
edges <- network_df %>%
    filter(source %in% kept_genes & target %in% kept_genes) %>%
    filter(PCC >= quantile(PCC, 0.25))

g <- igraph::graph_from_data_frame(
    edges %>% select(source, target, PCC),
    directed = FALSE,
    vertices = data.frame(name = kept_genes)
)

# Drop isolated nodes
isolated <- which(igraph::degree(g) == 0)
if (length(isolated) > 0) {
    g <- igraph::delete_vertices(g, isolated)
    cat("Removed", length(isolated), "isolated genes\n")
}

# Node attributes
gene_to_mod <- modules_df %>%
    filter(module_label %in% keep_modules) %>%
    select(name, module_label)
mod_lookup <- setNames(as.character(gene_to_mod$module_label), gene_to_mod$name)

V(g)$module <- mod_lookup[V(g)$name]
V(g)$domain <- domain_lookup[V(g)$module]
V(g)$degree <- igraph::degree(g)

# Edge attribute — must be set BEFORE create_layout
el <- igraph::ends(g, E(g))
is_within <- mod_lookup[el[,1]] == mod_lookup[el[,2]]
E(g)$weight <- ifelse(is_within, E(g)$PCC * 0.75, E(g)$PCC * 1.50)


# ==============================================================================
# Layout
# ==============================================================================
set.seed(428)
bb <- layout_as_backbone(g, keep = 0.1)

layout <- create_layout(g, layout = "stress")
layout$module <- V(g)$module
layout$domain <- V(g)$domain
layout$module_label <- label_lookup[layout$module]

# ==============================================================================
# Hull data (filtered to exclude outlier genes per module)
# ==============================================================================
hull_data <- as.data.frame(layout[, c("x", "y", "module_label")])
hull_data <- hull_data %>%
    group_by(module_label) %>%
    mutate(
        cx = mean(x),
        cy = mean(y),
        dist = sqrt((x - cx)^2 + (y - cy)^2),
        is_core = dist <= quantile(dist, 0.85)
    ) %>%
    ungroup() %>%
    filter(is_core)

# ==============================================================================
# Plot
# ==============================================================================
module_colors <- setNames(pal[domain_lookup[as.character(keep_modules)]],
                           label_lookup[as.character(keep_modules)])

# Vectorized within-module degree
kept_mods <- modules_df %>%
    filter(module_label %in% keep_modules, name %in% V(g)$name)

within_deg <- sapply(seq_len(nrow(kept_mods)), function(i) {
    gene <- kept_mods$name[i]
    mod_genes <- kept_mods$name[kept_mods$module_label == kept_mods$module_label[i]]
    sum((network_df$source == gene & network_df$target %in% mod_genes) |
        (network_df$target == gene & network_df$source %in% mod_genes))
})

kept_mods$within_degree <- within_deg

top_genes <- kept_mods %>%
    group_by(module_label) %>%
    slice_max(within_degree, n = 5, with_ties = FALSE) %>%
    ungroup()

top_gene_names <- top_genes$name

p <- ggraph(layout) +
    # Edges
    geom_edge_link0(
        aes(edge_alpha = is_within),
        edge_colour = "grey50",
        edge_width = 0.15
    ) +
    scale_edge_alpha_manual(values = c("FALSE" = 0.05, "TRUE" = 0.3), guide = "none") +
    # Hulls from filtered data
    geom_mark_hull(
        data = hull_data,
        aes(x = x, y = y, group = module_label, fill = module_label, label = module_label),
        concavity = 4,
        expand = unit(2, "mm"),
        alpha = 0.2,
        colour = NA,
        show.legend = FALSE,
        label.fontsize = 15,
    ) +
    # Nodes
    geom_node_point(
        aes(fill = domain, size = degree),
        shape = 21,
        colour = "grey30",
        stroke = 0.1,
        alpha = 0.85
    ) +
    # Gene labels (top 10% degree)
    geom_node_text(
        aes(filter = name %in% top_gene_names, label = name),
        size = 5,
        repel = TRUE,
        max.overlaps = 50,
        segment.color = "grey70",
        segment.size = 0.2,
        family = "sans"
    ) +
    # Colors
    scale_fill_manual(
        values = c(pal, module_colors),
        guide = guide_legend(title = "Domain", override.aes = list(size = 4, shape = 21))
    ) +
    scale_size(range = c(1, 6), guide = "none") +
    # Theme
    theme_graph(base_family = "sans") +
    theme(
        legend.position = "right",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey50")
    ) +
    labs(
        title = "Spatial Co-expression Network",
        subtitle = "10 domain-specific modules | PCC 0.8, Power 9"
    )

ggsave(file.path(plots_dir, "full_network_hulls_selected.pdf"),
       p, width = 16, height = 14, device = "pdf")
ggsave(file.path(plots_dir, "full_network_hulls_selected.png"),
       p, width = 16, height = 14, dpi = 300)

cat("Saved: full_network_hulls_selected.pdf / .png\n")
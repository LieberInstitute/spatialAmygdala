library(here)
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(ggplot2)

# ==============================================================================
# Configuration
# ==============================================================================
results_dir <- here("processed-data", "Visium", "14_smoothie_coexpression", "results")
plots_dir <- here("plots", "Visium", "14_smoothie_coexpression")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

network_df <- read.csv(file.path(results_dir, "network", "network_PCC0.9_clustPow9.csv"))
modules_df <- read.csv(file.path(results_dir, "modules_df_0.9_9.csv"))

MIN_GENES <- 5  # minimum genes per module to plot

# ==============================================================================
# Filter modules
# ==============================================================================
mod_sizes <- modules_df %>%
    group_by(module_label) %>%
    summarise(n_genes = n(), .groups = "drop") %>%
    filter(n_genes >= MIN_GENES) %>%
    arrange(desc(n_genes))

cat("Modules with >=", MIN_GENES, "genes:", nrow(mod_sizes), "\n")

# ==============================================================================
# Plot each module as a network graph
# ==============================================================================
pdf(file.path(plots_dir, "module_networks_0.9_9.pdf"), width = 10, height = 10)

for (i in seq_len(nrow(mod_sizes))) {
    mod <- mod_sizes$module_label[i]
    n_genes <- mod_sizes$n_genes[i]

    # Get genes in this module
    mod_genes <- modules_df %>%
        filter(module_label == mod) %>%
        pull(name)

    # Subset edges to within-module
    edges <- network_df %>%
        filter(source %in% mod_genes & target %in% mod_genes)

    if (nrow(edges) == 0) {
        cat("  Module", mod, "- no edges, skipping\n")
        next
    }

    # Build igraph object
    g <- graph_from_data_frame(
        edges %>% select(source, target, PCC),
        directed = FALSE,
        vertices = data.frame(name = mod_genes)
    )

    # Add degree as node attribute
    V(g)$degree <- degree(g)

    # Convert to tidygraph
    tg <- as_tbl_graph(g)

    # Choose layout based on size
    layout <- if (n_genes > 50) "fr" else if (n_genes > 20) "kk" else "stress"

    p <- ggraph(tg, layout = layout) +
        geom_edge_link(
            aes(alpha = PCC, width = PCC),
            color = "grey60",
            show.legend = FALSE
        ) +
        scale_edge_width(range = c(0.2, 1.5)) +
        scale_edge_alpha(range = c(0.15, 0.6)) +
        geom_node_point(
            aes(size = degree),
            color = "#2E75B6",
            alpha = 0.8
        ) +
        geom_node_text(
            aes(label = name),
            size = max(2, min(3.5, 80 / n_genes)),
            repel = TRUE,
            max.overlaps = 30,
            segment.color = "grey80",
            segment.size = 0.2
        ) +
        scale_size(range = c(2, 8), guide = "none") +
        labs(
            title = paste0("Module ", mod, "  (", n_genes, " genes, ", nrow(edges), " edges)"),
            subtitle = paste0("PCC cutoff = 0.8, power = 9")
        ) +
        theme_void() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey50"),
            plot.margin = margin(10, 10, 10, 10)
        )

    print(p)
    cat("  Module", mod, "(", n_genes, "genes,", nrow(edges), "edges )\n")
}

dev.off()
cat("\nSaved: module_networks_0.8_9.pdf\n")
library(here)
library(igraph)
library(ggraph)
library(graphlayouts)
library(dplyr)
library(patchwork)

# ==============================================================================
# Configuration
# ==============================================================================
results_dir <- here("processed-data", "Visium", "14_smoothie_coexpression", "results")
code_dir <- here("code", "Visium", "14_Smoothie_coexpression")
plots_dir <- here("plots", "Visium", "14_smoothie_coexpression")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

network_df <- read.csv(file.path(results_dir, "network", "network_PCC0.8_clustPow9.csv"))
modules_df <- read.csv(file.path(results_dir, "modules_df_0.8_9.csv"))

# Load marker genes
markers_df <- read.csv(here(code_dir, "top100_markers_bs_final_ITC_smoothed_wBLVM.csv"))
marker_genes <- unique(markers_df$gene)

# Load TF list
tf_df <- read.csv(file.path(code_dir, "human_TFs_Lambert2018.csv"))
tf_genes <- tf_df$HGNC.symbol[tf_df$Is.TF. == "Yes"]

cat("Marker genes:", length(marker_genes), "\n")
cat("TF genes:", length(tf_genes), "\n")

# ==============================================================================
# Modules to plot
# ==============================================================================
plot_modules <- c(
    "6"  = "LA",
    "27" = "BLD",
    "5"  = "PL",
    "10" = "CoA"
    #"13" = "CeA",
    #"19" = "MeA",
    #"20" = "AI",
    #"16" = "HPC"
)

# Domain colors for titles
pal <- c(
    AI  = "#D62728",
    BLD = "#9B59B6",
    PL  = "#f1e438ff",
    LA  = "#F4B400",
    CoA = "#5DA5DA",
    CeA = "#197d43ff",
    MeA = "#baf739ff",
    HPC = "#d6a8f8ff"
)

# ==============================================================================
# Build individual module network plots
# ==============================================================================
plot_list <- list()

for (mod in names(plot_modules)) {
    domain <- plot_modules[mod]
    mod_int <- as.integer(mod)

    # Get module genes
    mod_genes <- modules_df$name[modules_df$module_label == mod_int]

    # Subset edges: at least one end must be a TF
    edges <- network_df %>%
        filter(source %in% mod_genes & target %in% mod_genes) %>%
        filter(source %in% tf_genes | target %in% tf_genes)

    if (nrow(edges) == 0) {
        cat("Module", mod, "- no TF edges, skipping\n")
        next
    }

    # Only keep genes that appear in TF-connected edges
    edge_genes <- unique(c(edges$source, edges$target))

    g <- igraph::graph_from_data_frame(
        edges %>% select(source, target, PCC),
        directed = FALSE,
        vertices = data.frame(name = edge_genes)
    )

    # Node attributes
    V(g)$degree <- igraph::degree(g)
    domain_markers <- markers_df$gene[markers_df$cluster == domain]
    V(g)$is_marker <- V(g)$name %in% domain_markers
    V(g)$is_tf <- V(g)$name %in% tf_genes

    V(g)$gene_type <- "Target"
    V(g)$gene_type[V(g)$is_marker] <- "Marker"
    V(g)$gene_type[V(g)$is_tf] <- "TF"

    # Labels: all TFs, markers only if few enough
    n_markers <- sum(V(g)$is_marker)
    V(g)$show_label <- V(g)$is_tf | (V(g)$is_marker & n_markers <= 15)

    # Layout
    set.seed(42)
    layout <- create_layout(g, layout = "stress")

    p <- ggraph(layout) +
        geom_edge_link0(
            edge_colour = "grey80",
            edge_width = 0.4,
            edge_alpha = 0.4
        ) +
        # Target genes (non-TF, non-marker)
        geom_node_point(
            aes(filter = gene_type == "Target", size = degree),
            fill = "grey60",
            shape = 21,
            colour = "grey40",
            stroke = 0.5,
            alpha = 0.7
        ) +
        # Marker genes
        geom_node_point(
            aes(filter = gene_type == "Marker", size = degree),
            fill = "#1f4e79",
            shape = 21,
            colour = "grey30",
            stroke = 0.5,
            alpha = 0.9
        ) +
        # TF genes
        geom_node_point(
            aes(filter = gene_type == "TF", size = degree * 1.5),
            fill = "#c0392b",
            shape = 23,
            colour = "grey30",
            stroke = 0.4,
            alpha = 0.95
        ) +
        geom_node_text(
            aes(filter = gene_type == "TF", label = name),
            size = 5,
            colour = "black",
            fontface = "bold",
            repel = TRUE,
            max.overlaps = 30,
            segment.color = "grey50",
            segment.size = 0.2
        ) +
        scale_size(range = c(1.5, 8), guide = "none") +
        labs(title = paste0("M", mod, " - ", domain)) +
        theme_graph() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 22, face = "bold",
                                      colour = pal[domain]),
            plot.margin = margin(2, 2, 2, 2)
        )

    plot_list[[mod]] <- p
    cat("Module", mod, "-", domain, ":", igraph::vcount(g), "genes,",
        sum(V(g)$is_marker), "markers,", sum(V(g)$is_tf), "TFs\n")
}


# ==============================================================================
# Assemble patchwork (1 row x 8 panels)
# ==============================================================================
# 2 rows x 4 columns
combined <- wrap_plots(plot_list, nrow = 2)

ggsave(file.path(plots_dir, "module_networks_individual.pdf"),
       combined, width =  10, height = 10)
ggsave(file.path(plots_dir, "module_networks_individual.png"),
       combined, width = 10, height = 10, dpi = 300)

cat("\nSaved: module_networks_individual_row.pdf / .png\n")

# Print legend info
cat("\nLegend:\n")
cat("  Grey circles = other genes\n")
cat("  Dark blue circles = marker genes\n")
cat("  Red diamonds = transcription factors\n")
cat("  Labels shown for markers, TFs, and top 15% degree hub genes\n")
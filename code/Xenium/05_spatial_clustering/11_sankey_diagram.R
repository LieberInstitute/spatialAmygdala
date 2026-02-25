#!/usr/bin/env Rscript

# Sankey diagram of cluster curation iterations
# Install ggsankey if needed: devtools::install_github("davidsjoberg/ggsankey")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsankey)  # devtools::install_github("davidsjoberg/ggsankey")
  library(dplyr)
  library(tidyr)
})

output_dir <- here("plots", "Xenium", "05_spatial_clustering")

# ======== Define all iteration mappings ========
# Each list maps cluster_name -> vector of BANKSY res2 subclusters

V1 <- list(
  "C1"  = c(21,39,4,7),
  "C2"  = c(19,22,5,16),
  "C3"  = c(44,41,37,10,14),
  "C4"  = c(47),
  "C5"  = c(18,42,17,29,15,35),
  "C6"  = c(13,32,8,23,33,2,20,43),
  "C7"  = c(6,1,30),
  "C8"  = c(9,40,11,12,31),
  "C9"  = c(25),
  "C10" = c(45,46,3,36,34,26,27),
  "C11" = c(24,28,38),
  "C12" = c(52,48,49),
  "C13" = c(51,54),
  "C14" = c(53,50,55),
  "C15" = c(56,58,57,60),
  "C16" = c(64,65),
  "C17" = c(59,62,61,63)
)

V2 <- list(
  "C1"  = c(21,39,4,7),
  "C2"  = c(19,22,5,16),
  "C3"  = c(44,41,37,10,14),
  "C4"  = c(47),
  "C5"  = c(18,42,17,29,15,35),
  "C6"  = c(13,32,8,23,33,2,20,43),
  "C7"  = c(6,1,30),
  "C8"  = c(9),
  "C10" = c(11),
  "C11" = c(12,40),
  "C12" = c(31),
  "C13" = c(25),
  "C14" = c(45,46,3,36,34,26,27),
  "C15" = c(24,28,38),
  "C16" = c(52,48,49),
  "C17" = c(51,54),
  "C18" = c(53,50,55),
  "C19" = c(56,58,57,60),
  "C20" = c(64,65),
  "C21" = c(59,62,61,63)
)

V3 <- list(
  "C1"  = c(21,39,4,7),
  "C2"  = c(19),
  "C3"  = c(22),
  "C6"  = c(44,41,37,10,14),
  "C7"  = c(47),
  "C8"  = c(18,42,17,29,15,35),
  "C10" = c(32,8,23,33),
  "C11" = c(2,20,43),
  "C12" = c(6,1,30,5,16,13),
  "C13" = c(9),
  "C14" = c(11),
  "C15" = c(12,40),
  "C16" = c(31),
  "C17" = c(25),
  "C18" = c(45,46,3,36,34,26,27),
  "C19" = c(24,28,38),
  "C20" = c(52,48,49),
  "C21" = c(51,54),
  "C22" = c(53,50,55),
  "C23" = c(56,58,57,60),
  "C24" = c(64,65),
  "C25" = c(59,62,61,63)
)

V4 <- list(
  "C1"  = c(21,39,4,7),
  "C2"  = c(19),
  "C3"  = c(22),
  "C4"  = c(44,41,37,10,14),
  "C5"  = c(47),
  "C6"  = c(18,42,17,29,15,35),
  "C7"  = c(32),
  "C9"  = c(23),
  "C10" = c(33),
  "C11" = c(2,20,43,8),
  "C12" = c(6,1,30,5,16,13),
  "C13" = c(9),
  "C14" = c(11),
  "C15" = c(12,40),
  "C16" = c(31),
  "C17" = c(25),
  "C18" = c(45,46,3,36,34,26,27),
  "C19" = c(24,28,38),
  "C20" = c(52,48,49),
  "C21" = c(51,54),
  "C22" = c(53,50,55),
  "C23" = c(56,58,57,60),
  "C24" = c(64,65),
  "C25" = c(59,62,61,63)
)

V5 <- list(
  "C1"  = c(21,39,4,7),
  "C2"  = c(19),
  "C3"  = c(22),
  "C4"  = c(44,41,37,10,14),
  "C5"  = c(47),
  "C6"  = c(18,42,17,29,15,35),
  "C7"  = c(32),
  "C8"  = c(23),
  "C9"  = c(33),
  "C10" = c(2,20,43,8),
  "C11" = c(6,1,30,5,16,13),
  "C12" = c(9),
  "C13" = c(11),
  "C14" = c(12,40),
  "C15" = c(31),
  "C16" = c(25),
  "C17" = c(45),
  "C18" = c(46),
  "C19" = c(3),
  "C20" = c(36),
  "C21" = c(34),
  "C22" = c(26),
  "C23" = c(27),
  "C24" = c(24,28,38),
  "C25" = c(52,48,49),
  "C26" = c(51,54),
  "C27" = c(53,50,55),
  "C28" = c(56,58,57,60),
  "C29" = c(64,65),
  "C30" = c(59,62,61,63)
)

V6 <- list(
  "C1"  = c(21,39,4,7,26,19,22,44,41,37,10,14),
  "C2"  = c(47),
  "C3"  = c(18,42,17,29,15,35),
  "C4"  = c(32),
  "C5"  = c(23),
  "C6"  = c(33),
  "C7"  = c(2,20,43,8),
  "C8"  = c(6,1,30,5,16,13),
  "C9"  = c(9),
  "C10" = c(11),
  "C11" = c(12,40),
  "C12" = c(31),
  "C13" = c(25),
  "C14" = c(45),
  "C15" = c(46),
  "C16" = c(3),
  "C17" = c(36),
  "C18" = c(34,27),
  "C19" = c(24,28,38),
  "C20" = c(52,48,49),
  "C21" = c(51,54),
  "C22" = c(53,50,55),
  "C23" = c(56,58,57,60),
  "C24" = c(64,65),
  "C25" = c(59,62,61,63)
)

# ======== Build long-form data ========
# Each row = one BANKSY subcluster, with its cluster assignment at each iteration

all_iterations <- list(V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5, V6 = V6)

# Reverse lookup: subcluster -> cluster name, per iteration
reverse_lookup <- function(mapping) {
  out <- list()
  for (cl_name in names(mapping)) {
    for (sub in mapping[[cl_name]]) {
      out[[as.character(sub)]] <- cl_name
    }
  }
  out
}

rev_maps <- lapply(all_iterations, reverse_lookup)
all_subs <- sort(unique(as.integer(unlist(lapply(all_iterations, function(x) unlist(x))))))

# Build wide data frame: one row per subcluster
df_wide <- data.frame(subcluster = all_subs)
for (iter_name in names(rev_maps)) {
  df_wide[[iter_name]] <- sapply(as.character(df_wide$subcluster), function(s) {
    rev_maps[[iter_name]][[s]] %||% NA_character_
  })
}

# Helper: extract numeric part from cluster name for sorting
cluster_num <- function(x) as.integer(gsub(".*C(\\d+)$", "\\1", x))

# Prefix each column so node names are unique across stages,
# and build numerically-sorted factor levels per iteration
for (iter_name in names(all_iterations)) {
  df_wide[[iter_name]] <- paste0(iter_name, ": ", df_wide[[iter_name]])
  # Get unique values, sort by the numeric cluster number
  uvals <- unique(df_wide[[iter_name]])
  uvals <- uvals[order(cluster_num(uvals))]
  df_wide[[iter_name]] <- factor(df_wide[[iter_name]], levels = uvals)
}

# Convert to ggsankey long format
df_long <- df_wide %>%
  make_long(V1, V2, V3, V4, V5, V6)

# Build a single global factor level vector: all iterations, each sorted numerically
all_node_levels <- unlist(lapply(names(all_iterations), function(iter_name) {
  uvals <- unique(grep(paste0("^", iter_name, ":"), df_long$node, value = TRUE))
  uvals[order(cluster_num(uvals))]
}))

# Apply to both node and next_node so ggsankey respects ordering
df_long$node      <- factor(df_long$node,      levels = all_node_levels)
df_long$next_node <- factor(df_long$next_node,  levels = all_node_levels)

# Strip prefix for display labels
df_long$label <- gsub("^V[0-9]+: ", "", as.character(df_long$node))

# ======== Detect changed nodes ========
# A node is "changed" if its subcluster composition differs from the
# same-named cluster in the previous OR next iteration.
# Compare by sorting subclusters so order doesn't matter.

iter_names <- names(all_iterations)

# For each iteration, build a named list: cluster_name -> sorted subcluster signature
iter_signatures <- lapply(all_iterations, function(mapping) {
  lapply(mapping, function(subs) paste(sort(subs), collapse = ","))
})

# For each node (iter:cluster), check if it changed vs neighbors
changed_nodes <- character(0)

for (i in seq_along(iter_names)) {
  iter <- iter_names[i]
  sigs <- iter_signatures[[i]]
  
  for (cl_name in names(sigs)) {
    node_id <- paste0(iter, ": ", cl_name)
    this_sig <- sigs[[cl_name]]
    is_changed <- FALSE
    
    # Compare to previous iteration: does any cluster in prev have the same signature?
    if (i > 1) {
      prev_sigs <- unlist(iter_signatures[[i - 1]])
      if (!this_sig %in% prev_sigs) {
        is_changed <- TRUE
      }
    }
    
    # Compare to next iteration: does any cluster in next have the same signature?
    if (i < length(iter_names)) {
      next_sigs <- unlist(iter_signatures[[i + 1]])
      if (!this_sig %in% next_sigs) {
        is_changed <- TRUE
      }
    }
    
    # V1 is the starting point — mark all as changed for initial context
    if (i == 1) is_changed <- TRUE
    
    if (is_changed) changed_nodes <- c(changed_nodes, node_id)
  }
}

# Assign fill group: changed nodes get their own color, unchanged get "unchanged"
df_long$fill_group <- ifelse(as.character(df_long$node) %in% changed_nodes,
                              as.character(df_long$node),
                              "unchanged")

# Build color palette: distinct hues for changed, grey for unchanged
changed_unique <- unique(df_long$fill_group[df_long$fill_group != "unchanged"])
n_changed <- length(changed_unique)
changed_colors <- setNames(
  scales::hue_pal()(n_changed),
  changed_unique
)
all_colors <- c(changed_colors, "unchanged" = "grey80")

# ======== Plot ========
p <- ggplot(df_long, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                          fill = fill_group, label = label)) +
  geom_sankey(flow.alpha = 0.35, node.color = "grey30", show.legend = FALSE) +
  geom_sankey_label(size = 2.5, color = "white", fill = "grey20", alpha = 0.7) +
  scale_fill_manual(values = all_colors) +
  scale_x_discrete(labels = c("V1\nInitial\n(17 cl.)",
                               "V2\nSplit C7\n(20 cl.)",
                               "V3\nUncollapse C2\n(22 cl.)",
                               "V4\nUncollapse C10\n(24 cl.)",
                               "V5\nUncollapse C18\n(30 cl.)",
                               "V6\nFinal collapse\n(25 cl.)")) +
  theme_sankey(base_size = 12) +
  labs(title = "Cluster Curation: Sankey of Manual Merging/Splitting",
       subtitle = "Colored = changed composition vs neighboring iteration | Grey = unchanged",
       x = NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10),
    axis.text.x = element_text(size = 9, lineheight = 1.1)
  )

ggsave(here(output_dir, "cluster_curation_sankey.pdf"), p, width = 18, height = 14)

suppressPackageStartupMessages({
    library("here")
    library("ComplexHeatmap")
    library("circlize")
    library("grid")
})

# ---------------------------
# Directories  (mirrors config.sh: WORKDIR, DOWNSTREAM_DIR, FIG_DIR, ANNOT_COL)
# ---------------------------
workdir   <- here("processed-data", "Visium", "16_LDSC", "scDRS")
down_dir  <- file.path(workdir, "downstream")   # SUB=pilot -> "downstream/pilot"
plots_dir <- file.path(workdir, "figures")      # $FIG_DIR
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

annot <- "BS_k16_Semisupervised_wAI"   # $ANNOT_COL

# ---------------------------
# Traits to keep, and their display names
# ---------------------------
trait_labels <- c(
    GSCAN_SmkInit   = "Smoking initiation",
    Neuroticism     = "Neuroticism",
    BMI             = "BMI",
    EduYears        = "Education years",
    Anorexia        = "Anorexia",
    SCZ             = "Schizophrenia",
    ADHD            = "ADHD",
    Insomnia        = "Insomnia",
    Intelligence    = "Intelligence",
    MDD             = "Depression",
    GSCAN_DrnkWk    = "Drinks per week",
    BIP_2024        = "Bipolar disorder",
    Alzheimer_V3    = "Alzheimer's",
    Autism          = "Autism",
    PD              = "Parkinson's",
    PTSD_F3         = "PTSD",
    Stroke_2022_Any = "Stroke",
    OUD_META        = "Opioid use disorder",
    Epilepsy_GGE    = "Epilepsy (GGE)",
    Height          = "Height",
    T2D             = "Type 2 diabetes"
)

# ---------------------------
# Color palette
# ---------------------------
pal <- c(
    IA          = "#D62728",
    BM          = "#E67E22",
    BLD         = "#9B59B6",
    PL          = "#f1e438",
    BL          = "#035185",
    LA          = "#F4B400",
    CoA         = "#5DA5DA",
    CeA         = "#197d43",
    MeA         = "#baf739",
    HPC         = "#d6a8f8",
    CHAT        = "#A0522D",
    Endothelial = "#444444",
    WM.1        = "#BBBBBB",
    WM.2        = "#DDDDDD",
    CLA         = "#FF69B4"
)

# ---------------------------
# Domain order and exclusions
# ---------------------------
domain_order <- c(
    "CLA", "HPC", "BLD", "LA", "PL", "BL", "BM", "CoA",
    "IA", "CeA", "MeA", "CHAT", "Endothelial", "WM.1"
)

domain_drop <- "WM.2"

# ---------------------------
# Load group analysis
# ---------------------------
files <- sort(Sys.glob(file.path(down_dir, paste0("*.scdrs_group.", annot))))
length(files)

grp <- do.call(rbind, lapply(files, function(f) {
    d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
    names(d)[1] <- "domain"
    d$trait <- sub(paste0("\\.scdrs_group\\.", annot, "$"), "", basename(f))
    d
}))

# strip the " (n/7)" donor suffix if it's baked into the labels
grp$domain <- trimws(sub("\\s*\\(\\d+/\\d+\\)\\s*$", "", grp$domain))

# AI -> IA (intercalated)
grp$domain[grp$domain == "AI"] <- "IA"

head(grp)
sort(unique(grp$trait))

# ---------------------------
# Subset to the traits we want
# ---------------------------
missing <- setdiff(names(trait_labels), unique(grp$trait))
if (length(missing)) message("no group file for: ", paste(missing, collapse = ", "))

grp <- grp[grp$trait %in% names(trait_labels), ]
grp$trait <- trait_labels[grp$trait]

# ---------------------------
# Stars (BH FDR < 0.05 across the plotted grid)
# ---------------------------
grp <- grp[!grp$domain %in% domain_drop, ]

grp$fdr  <- p.adjust(grp$assoc_mcp, method = "BH")
grp$star <- ifelse(grp$fdr < 0.05, "*", "")

table(grp$star)

# ---------------------------
# Long -> matrix; fixed row order, columns by mean score
# ---------------------------
z_mat <- tapply(grp$assoc_mcz, list(grp$domain, grp$trait), identity)
s_mat <- tapply(grp$star,      list(grp$domain, grp$trait), identity)
s_mat[is.na(s_mat)] <- ""

extra <- setdiff(rownames(z_mat), domain_order)
if (length(extra)) message("domain not in domain_order, dropped: ",
                           paste(extra, collapse = ", "))

domain_order <- intersect(domain_order, rownames(z_mat))
trait_order  <- names(sort(colMeans(z_mat, na.rm = TRUE), decreasing = TRUE))

z_mat <- z_mat[domain_order, trait_order]
s_mat <- s_mat[domain_order, trait_order]

# ---------------------------
# Annotation and color scale
# ---------------------------
dom_cols <- pal[domain_order]

row_ha <- rowAnnotation(
    Domain = factor(domain_order, levels = domain_order),
    col = list(Domain = dom_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(4, "mm")
)

# n_cell is identical across traits for a given domain, so take the first.
# c() strips the 1-D array dim tapply leaves behind -- anno_barplot reads a
# non-null dim as "matrix" and errors in rowSums().
n_spots <- c(tapply(grp$n_cell, grp$domain, function(x) x[1]))[domain_order]
n_spots

row_bar <- rowAnnotation(
    `Spots` = anno_barplot(
        n_spots,
        gp = gpar(fill = dom_cols, col = NA),
        border = FALSE,
        axis_param = list(gp = gpar(fontsize = 8)),
        width = unit(2.2, "cm")
    ),
    annotation_name_gp = gpar(fontsize = 10),
    annotation_name_rot = 0,
    annotation_name_side = "top"
)

z_lim  <- max(abs(z_mat), na.rm = TRUE)
z_cols <- colorRamp2(c(-z_lim, 0, z_lim), c("#2166AC", "white", "#B2182B"))

# ---------------------------
# Heatmap
# ---------------------------
ht <- Heatmap(
    z_mat,
    name = "Disease relevance",
    col = z_cols,
    left_annotation = row_ha,
    right_annotation = row_bar,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    rect_gp = gpar(col = "white", lwd = 0.5),
    na_col = "grey92",
    cell_fun = function(j, i, x, y, width, height, fill) {
        if (nzchar(s_mat[i, j])) {
            grid::grid.text("*", x, y, vjust = 0.75,
                            gp = grid::gpar(fontsize = 11, col = "white"))
        }
    },
    width  = unit(7 * ncol(z_mat), "mm"),
    height = unit(7 * nrow(z_mat), "mm"),
    heatmap_legend_param = list(
        title = "Disease\nrelevance (z)",
        legend_height = unit(3.5, "cm")
    )
)

# ---------------------------
# Draw
# ---------------------------
w <- 5.8 + 0.29 * ncol(z_mat)
h <- 2.5 + 0.29 * nrow(z_mat)

pdf(file.path(plots_dir, "scdrs_domain_heatmap.pdf"), width = w, height = h)
draw(ht,
     column_title = "scDRS domain-level association (* BH FDR < 0.05)",
     column_title_gp = gpar(fontsize = 13),
     heatmap_legend_side = "right")
dev.off()

write.table(grp, file.path(plots_dir, "scdrs_domain_trait_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Saved to: ", plots_dir)
library("SpatialFeatureExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("ggspavis")
library("scCustomize")
library("ggplot2")
library("scattermore")

# ======= Load base spe =======
spe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control", "spe_qc_016.rds"))

# drop out of tissue for all except Br8325_MeA
spe <- spe[, spe$in_tissue | spe$sample_id == "Br8325_MeA"]

# ======= Load Banksy output for resolution 0.6 =======
banksy_dir <- here("processed-data", "VisiumHD", "04_clustering", "016_wAI_markers")
banksy_file <- file.path(banksy_dir, "Banksy_integrated_lambda_0.8_res0.6.rds")

if (!file.exists(banksy_file)) {
    stop("Banksy file for resolution 0.6 not found: ", banksy_file)
}

cat("Loading Banksy output for resolution 0.6\n")
sfe.sub <- readRDS(banksy_file)

# Find the cluster column
clust_cols <- grep("clust_HARMONY_M1_lam0\\.8", colnames(colData(sfe.sub)), value = TRUE)
if (length(clust_cols) == 0) {
    stop("No cluster column found matching pattern clust_HARMONY_M1_lam0.8")
}
clust_col <- clust_cols[1]
cat("Cluster column:", clust_col, "\n")

# Transfer labels to full spe
spe[[clust_col]] <- factor(sfe.sub[[clust_col]])

# ======= Define cluster renaming =======
default_labels <- c(
    "1"  = "MeA",
    "2"  = "BM",
    "3"  = "Neuropil",
    "4"  = "BLD",
    "5"  = "BM",
    "6"  = "WM",
    "7"  = "Ventricle",
    "8"  = "BLD",
    "9"  = "CeA",
    "10" = "CoA",
    "11" = "Neuropil",
    "12" = "BM",
    "13" = "IA",
    "14" = "Vascular",
    "15" = "DELETE",
    "16" = "BM",
    "17" = "BM",
    "18" = "Neuropil",
    "19" = "CHAT",
    "20" = "DELETE",
    "21" = "DELETE",
    "22" = "DELETE"
)

sample_overrides <- list(
    "Br9280_ITC" = c("2" = "LA"),
    "Br8325_CeA" = c("12" = "BLD")
)

# fail loudly rather than silently skipping an override that matches nothing
.ids <- unique(as.character(spe$sample_id))
.bad <- setdiff(names(sample_overrides), .ids)
if (length(.bad) > 0) {
    stop("sample_overrides target sample_id(s) not present: ", paste(.bad, collapse = ", "),
         "\n  available: ", paste(.ids, collapse = ", "))
}
for (.sid in names(sample_overrides)) {
    .present <- unique(as.character(spe[[clust_col]])[spe$sample_id == .sid])
    .miss <- setdiff(names(sample_overrides[[.sid]]), .present)
    if (length(.miss) > 0) {
        warning("cluster(s) ", paste(.miss, collapse = ", "), " not found in ", .sid,
                " -- override will not apply")
    }
}

# ======= Apply renaming =======
original_clusters <- as.character(spe[[clust_col]])
sample_ids <- as.character(spe$sample_id)
new_labels <- character(length(original_clusters))

for (i in seq_along(original_clusters)) {
    cl <- original_clusters[i]
    sid <- sample_ids[i]
    if (sid %in% names(sample_overrides) && cl %in% names(sample_overrides[[sid]])) {
        new_labels[i] <- sample_overrides[[sid]][cl]
    } else {
        new_labels[i] <- default_labels[cl]
    }
}

# ======= Remove DELETE spots =======
keep <- new_labels != "DELETE"
cat("Removing", sum(!keep), "spots marked DELETE\n")
spe <- spe[, keep]
new_labels <- new_labels[keep]

label_col <- "spatial_domain"
spe[[label_col]] <- factor(new_labels)

cat("Final spatial domains:\n")
print(table(spe[[label_col]]))

# ======= Reorient specific samples =======
# FALSE renders the raw orientation, to a separate PDF, so you can pick values.
apply_transforms <- FALSE

# TRUE renders all 8 rigid orientations per sample to a contact sheet, so the
# correct one can be read straight off the panel label. Requires raw coords.
orientation_scout <- TRUE

if (orientation_scout && apply_transforms) {
    stop("set apply_transforms <- FALSE when orientation_scout is TRUE, ",
         "otherwise the scout starts from already-transformed coordinates")
}

# Applied in the order stated: rotate first, then mirror.
#   rotate: degrees COUNTER-clockwise as drawn on the page
#   mirror: "vertical"   = top-bottom flip (y -> -y)
#           "horizontal" = left-right flip (x -> -x)
#           "none"
sample_transforms <- list(
    "Br8325_CeA" = list(rotate =   0, mirror = "vertical"),
    "Br8325_MeA" = list(rotate =   0, mirror = "vertical"),
    "Br9280_CeA" = list(rotate = 270, mirror = "vertical"),
    "Br9280_ITC" = list(rotate = 270, mirror = "vertical"),
    "Br9280_MeA" = list(rotate = 270, mirror = "vertical")
)

.bad <- setdiff(names(sample_transforms), unique(as.character(spe$sample_id)))
if (length(.bad) > 0) {
    stop("sample_transforms target sample_id(s) not present: ", paste(.bad, collapse = ", "))
}

transform_xy <- function(xy, rotate = 0, mirror = c("none", "vertical", "horizontal")) {
    mirror <- match.arg(mirror)
    th  <- rotate * pi / 180
    p   <- sweep(xy, 2, colMeans(xy))
    out <- cbind(p[, 1] * cos(th) - p[, 2] * sin(th),
                 p[, 1] * sin(th) + p[, 2] * cos(th))
    if (mirror == "vertical")   out[, 2] <- -out[, 2]
    if (mirror == "horizontal") out[, 1] <- -out[, 1]
    out
}

if (apply_transforms) {
    for (s in names(sample_transforms)) {
        tf  <- sample_transforms[[s]]
        idx <- which(spe$sample_id == s)
        new <- transform_xy(spatialCoords(spe)[idx, 1:2, drop = FALSE],
                            rotate = tf$rotate, mirror = tf$mirror)
        spatialCoords(spe)[idx, 1] <- new[, 1]
        spatialCoords(spe)[idx, 2] <- new[, 2]
        cat("Reoriented", s, ": rotate", tf$rotate, "deg, mirror", tf$mirror,
            "(", length(idx), "spots )\n")
    }
} else {
    cat("apply_transforms = FALSE -- plotting raw orientation\n")
}

# ======= Color palette =======
domain_colors <- c(
    IA = "#D62728", BM = "#E67E22", BLD = "#9B59B6",
    LA = "#F4B400", CoA = "#5DA5DA", CeA = "#197d43ff",
    MeA = "#baf739ff", CHAT = "#A0522D",
    Ventricle = "#666666", Vascular = "#333333",
    WM = "#BBBBBB", Neuropil = "#DDDDDD"
)

# ======= Plot one sample per page using scattermore =======
plots_dir <- here("plots", "VisiumHD", "04_clustering", "016_wAI_markers")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

unique_samples <- sort(unique(spe$sample_id))

# ======= Orientation contact sheet =======
# One page per sample, all 8 rigid orientations, each panel titled with the
# exact rotate/mirror pair that produced it. Pick the right panel, then copy
# its label into sample_transforms above.
if (orientation_scout) {
    scout_grid <- expand.grid(
        rotate = c(0, 90, 180, 270),
        mirror = c("none", "vertical"),
        stringsAsFactors = FALSE
    )

    scout_path <- file.path(plots_dir, "orientation_contact_sheet.pdf")
    pdf(file = scout_path, width = 16, height = 8)

    for (s in names(sample_transforms)) {
        idx <- which(spe$sample_id == s)
        xy0 <- spatialCoords(spe)[idx, 1:2, drop = FALSE]
        dom <- droplevels(spe[[label_col]][idx])

        panels <- lapply(seq_len(nrow(scout_grid)), function(i) {
            g  <- scout_grid[i, ]
            xy <- transform_xy(xy0, rotate = g$rotate, mirror = g$mirror)
            ggplot(data.frame(x = xy[, 1], y = xy[, 2], domain = dom),
                   aes(x, y, color = domain)) +
                geom_scattermore(pointsize = 1.5, pixels = c(600, 600)) +
                scale_color_manual(values = domain_colors, na.value = "black") +
                coord_equal() +
                ggtitle(paste0("rotate = ", g$rotate, ", mirror = \"", g$mirror, "\"")) +
                theme_minimal(base_size = 10) +
                theme(legend.position = "none",
                      plot.title = element_text(hjust = 0.5, size = 11),
                      axis.title = element_blank(), axis.text = element_blank(),
                      panel.grid = element_blank(),
                      panel.border = element_rect(colour = "grey70", fill = NA))
        })

        print(wrap_plots(panels, nrow = 2) +
              plot_annotation(title = s,
                              theme = theme(plot.title = element_text(hjust = 0.5,
                                                                      size = 18,
                                                                      face = "bold"))))
        cat("Scouted:", s, "\n")
    }

    dev.off()
    cat("Contact sheet written to:", scout_path, "\n")
}

pdf_path <- file.path(plots_dir, if (apply_transforms) {
    "Banksy_integrated_lambda_0.8_res0.6_renamed_v4.pdf"
} else {
    "Banksy_integrated_lambda_0.8_res0.6_renamed_UNROTATED.pdf"
})
pdf(file = pdf_path, width = 10, height = 10)

for (s in unique_samples) {
    idx <- which(spe$sample_id == s)
    coords <- as.data.frame(spatialCoords(spe)[idx, ])
    colnames(coords) <- c("x", "y")
    # droplevels + drop = TRUE => legend lists only domains present in this sample
    coords$domain <- droplevels(spe[[label_col]][idx])

    p <- ggplot(coords, aes(x = x, y = y, color = domain)) +
        geom_scattermore(pointsize = 2, pixels = c(1024, 1024)) +
        scale_color_manual(values = domain_colors, drop = TRUE,
                           na.value = "black") +
        coord_equal() +
        ggtitle(s) +
        theme_minimal(base_size = 14) +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            legend.text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank()
        ) +
        guides(color = guide_legend(override.aes = list(size = 4)))

    print(p)
    cat("Plotted:", s, "\n")
}

dev.off()
cat("\nDone! Plot saved to:", pdf_path, "\n")



# save -- only when rotations are applied, so a scouting run can't overwrite
# the saved object with un-rotated coordinates
if (apply_transforms) {
    saveRDS(spe, here("processed-data", "VisiumHD", "04_clustering", "spe_banksy_renamed.rds"))
} else {
    cat("apply_transforms = FALSE -- skipping saveRDS\n")
}
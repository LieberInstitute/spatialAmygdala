library(jsonlite)
library(ggplot2)
library(ggrepel)

proj <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
data_dir <- file.path(proj, "code/Visium/16_LDSC/gsMap/figdata")
plots_dir <- file.path(proj, "processed-data/Visium/16_LDSC/gsMap/figures")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

P <- read.csv(file.path(data_dir, "test2_cauchy_matrix.csv"), row.names = 1, check.names = FALSE)
grp_df <- read.csv(file.path(data_dir, "trait_groups.csv"), row.names = 1)
grp <- setNames(grp_df$group, rownames(grp_df))
gsel <- read.csv(file.path(data_dir, "gwas_selection_table.csv"))

read_json <- function(f){
    x <- paste(readLines(f, warn = FALSE), collapse = "")
    x <- gsub("-Infinity|Infinity|NaN", "null", x)
    fromJSON(x)
}

pw_full <- read_json(file.path(data_dir, "chi2_full.json"))
pw_single <- read_json(file.path(data_dir, "chi2_single.json"))

picked <- gsel$trait[as.logical(gsel$selected)]
singles <- c("ADHD", "Autism", "Anorexia", "PTSD_F3", "Insomnia", "Neuroticism",
             "EduYears", "Intelligence", "PD", "GSCAN_AgeSmk", "GSCAN_CigDay",
             "GSCAN_SmkCes", "Height", "BMI", "T2D")
traits <- c(picked, singles[!singles %in% picked])
traits <- traits[traits %in% rownames(P)]

pow <- rbind(pw_full[pw_full$trait %in% picked, ], pw_single)
rownames(pow) <- pow$trait
stopifnot(all(traits %in% rownames(pow)))

p_mat <- as.matrix(P[traits, , drop = FALSE])
p_mat[p_mat < 1e-300] <- 1e-300

df <- data.frame(
    trait = traits,
    mean_chi2 = pow[traits, "mean_chi2"],
    mean_neglog10p = rowMeans(-log10(p_mat), na.rm = TRUE),
    group = grp[traits]
)

cols <- c(
    "Psychiatric" = "#3B6FB6",
    "Substance use" = "#4E9F50",
    "Cognitive/behavioral" = "#8E6FBF",
    "Neurological" = "#D97C2B",
    "Non-brain control" = "#6E6E6E"
)

# change PTSD_F3 to PTSD for display purposes
df$trait[df$trait == "PTSD_F3"] <- "PTSD"

highlight <- c("SCZ", "MDD", "PTSD", "Height")
ct <- cor.test(df$mean_chi2, df$mean_neglog10p, method = "spearman", exact = FALSE)
thr <- -log10(0.05)
hl <- df[df$trait %in% highlight, ]

p <- ggplot(df, aes(mean_chi2, mean_neglog10p)) +
    geom_hline(yintercept = thr, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_point(aes(color = group), size = 2.3, alpha = 0.5) +
    geom_point(data = hl, aes(fill = group), shape = 21, color = "white", stroke = 1, size = 4) +
    geom_text_repel(data = hl[hl$trait != "PTSD", ],
                    aes(label = trait, color = group), fontface = "bold",
                    size = 3.5, segment.color = "#888888", show.legend = FALSE) +
    geom_text_repel(data = hl[hl$trait == "PTSD", ],
                    aes(label = trait, color = group), fontface = "bold",
                    size = 3.5, nudge_y = 0.55, direction = "y",
                    segment.color = "#888888", show.legend = FALSE) +
    scale_x_log10(breaks = c(1, 2, 4, 8)) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols, guide = "none") +
    labs(x = expression("GWAS power (mean genome-wide "*chi^2*")"),
         y = expression("Mean "*-log[10](P)*" across domains"),
         title = "Enrichment magnitude tracks GWAS power",
         subtitle = paste0("Spearman \u03c1 = ", sprintf("%.2f", unname(ct$estimate)),
                           ", P = ", format(ct$p.value, scientific = TRUE, digits = 1),
                           ", n = ", nrow(df)), color = NULL) +
    guides(color = guide_legend(ncol = 1, override.aes = list(alpha = 1, size = 2.5))) +
    theme_classic(base_size = 11) +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.8), color = NA),
        legend.key.height = unit(0.4, "cm"),
        legend.spacing.y = unit(0.05, "cm")
    )

ggsave(file.path(plots_dir, "gsmap_power_vs_enrichment.pdf"), p, width = 3.5, height = 4)
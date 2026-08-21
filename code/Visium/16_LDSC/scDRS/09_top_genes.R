library("here")
library("dplyr")
library("ggplot2")
library("patchwork")

# Save directories
workdir <- here("processed-data", "Visium", "16_LDSC", "scDRS")
down_dir <- file.path(workdir, "downstream")
plots_dir <- file.path(workdir, "figures")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

annot <- "BS_k16_Semisupervised_wAI"

# Domain colors
pal <- c(
    IA = "#D62728",
    BM = "#E67E22",
    BLD = "#9B59B6",
    PL = "#f1e438",
    BL = "#035185",
    LA = "#F4B400",
    CoA = "#5DA5DA",
    CeA = "#197d43",
    MeA = "#baf739",
    HPC = "#d6a8f8",
    CHAT = "#A0522D",
    Endothelial = "#444444",
    WM.1 = "#BBBBBB",
    WM.2 = "#DDDDDD",
    CLA = "#FF69B4"
)

# ===== Load scDRS domain results =====

files <- sort(Sys.glob(
    file.path(down_dir, paste0("*.scdrs_group.", annot))
))

grp <- do.call(rbind, lapply(files, function(f) {
    d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
    names(d)[1] <- "domain"
    d$trait <- sub(paste0("\\.scdrs_group\\.", annot, "$"), "", basename(f))
    d
}))

grp$domain <- trimws(sub("\\s*\\(\\d+/\\d+\\)\\s*$", "", grp$domain))
grp$domain[grp$domain == "AI"] <- "IA"

# ===== H1: psychiatric enrichment by domain =====

psych_traits <- c(
    "SCZ",
    "BIP_2024",
    "MDD",
    "PTSD_F3",
    "ADHD",
    "Autism",
    "Anorexia",
    "OUD_META"
)

gm_domains <- c(
    "CLA", "HPC", "BLD", "LA", "PL", "BL",
    "BM", "CoA", "IA", "CeA", "MeA", "CHAT"
)

psych <- grp |>
    filter(
        trait %in% psych_traits,
        domain %in% gm_domains
    ) |>
    group_by(trait) |>
    mutate(
        enrichment_rank = percent_rank(assoc_mcz)
    ) |>
    ungroup()

# Average within-trait rank across psychiatric traits
psych_domains <- psych |>
    group_by(domain) |>
    summarize(
        mean_rank = mean(enrichment_rank, na.rm = TRUE),
        mean_z = mean(assoc_mcz, na.rm = TRUE),
        n_traits = n(),
        .groups = "drop"
    ) |>
    arrange(desc(mean_rank))

psych_domains

# order by consensus psychiatric enrichment
psych_domains$domain <- factor(
    psych_domains$domain,
    levels = rev(psych_domains$domain)
)

p1 <- ggplot(
    psych_domains,
    aes(x = mean_rank, y = domain)
) +
    geom_segment(
        aes(
            x = 0,
            xend = mean_rank,
            yend = domain,
            color = domain
        ),
        linewidth = 0.8
    ) +
    geom_point(
        aes(color = domain),
        size = 3
    ) +
    scale_color_manual(values = pal, guide = "none") +
    scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25)
    ) +
    labs(
        x = "Mean psychiatric enrichment rank",
        y = NULL,
        title = "Psychiatric enrichment by domain"
    ) +
    theme_classic(base_size = 10)

pdf(file.path(plots_dir, "scdrs_panel_H1.pdf"), width = 2.5,height = 4.5)
p1
dev.off()


 library("ggrepel")



# ============================================================
# H2: psychiatric enrichment across donors 
# ===========================================================



library("here")
library("dplyr")
library("ggplot2")

plots_dir <- here(
    "processed-data", "Visium", "16_LDSC", "scDRS", "figures"
)

pal <- c(
    IA = "#D62728",
    BM = "#E67E22",
    BLD = "#9B59B6",
    PL = "#f1e438",
    BL = "#035185",
    LA = "#F4B400",
    CoA = "#5DA5DA",
    CeA = "#197d43",
    MeA = "#baf739",
    HPC = "#d6a8f8",
    CHAT = "#A0522D",
    CLA = "#FF69B4"
)

# ===== Load donor-level enrichment =====

psych_donor <- read.delim(
    file.path(
        plots_dir,
        "psych_enrichment_donor_domain.tsv"
    )
)

psych_summary <- read.delim(
    file.path(
        plots_dir,
        "psych_enrichment_domain_summary.tsv"
    )
)

psych_summary

# order domains by mean across donors
domain_order <- psych_summary |>
    arrange(mean_psych_rank) |>
    pull(domain)

psych_donor$domain <- factor(
    psych_donor$domain,
    levels = domain_order
)

# show number of donors represented for each domain
domain_labels <- setNames(
    paste0(
        psych_summary$domain,
        "  (n=", psych_summary$n_donors, ")"
    ),
    psych_summary$domain
)

# ===== Plot =====

p <- ggplot(
    psych_donor,
    aes(
        x = mean_psych_rank,
        y = domain,
        color = domain
    )
) +
    geom_boxplot(
        width = 0.5,
        outlier.shape = NA,
        linewidth = 0.6
    ) +
    geom_jitter(
        height = 0.08,
        width = 0,
        size = 1.5,
        alpha = 0.5
    ) +
    scale_color_manual(
        values = pal,
        guide = "none"
    ) +
    scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25)
    ) +
    labs(
        x = "Mean rank",
        y = NULL,
        title = "Psychiatric enrichment"
    ) +
    theme_classic(base_size = 11)


pdf(
    file.path(
        plots_dir,
        "scdrs_psych_enrichment_across_donors.pdf"
    ),
    width = 2.5,
    height = 4
)

p

dev.off()





# ============================================================
# H3: psychiatric enrichment across donors 
# ===========================================================



library("here")
library("dplyr")
library("ggplot2")
library("ggrepel")

plots_dir <- here(
    "processed-data", "Visium", "16_LDSC", "scDRS", "figures"
)

# ===== Load scDRS contributions =====

contrib <- read.delim(
    file.path(
        plots_dir,
        "SCZ_LA_gene_contributions.tsv"
    )
)

# ===== Load spatialLIBD marker statistics =====

markers <- readRDS(
    here(
        "processed-data", "Visium", "08_marker_genes",
        "spatialLIBD_model_results.rds"
    )
)

enrichment <- as.data.frame(markers$enrichment)

if (!"gene" %in% colnames(enrichment)) {
    enrichment$gene <- rownames(enrichment)
}

la <- enrichment |>
    select(
        gene,
        logFC_LA
    )

# ===== Join =====

plot_dat <- contrib |>
    inner_join(la, by = "gene")

# top 10 genes by actual contribution to SCZ score in LA
top_genes <- plot_dat |>
    slice_max(
        contribution_domain,
        n = 10
    )

# ===== Plot =====

p <- ggplot(
    plot_dat,
    aes(
        x = magma_weight,
        y = logFC_LA,
        color = contribution_domain
    )
) +
    geom_point(
        size = 2,
        alpha = 0.7
    ) +
    geom_text_repel(
        data = top_genes,
        aes(label = gene),
        color = "black",
        fontface = "italic",
        size = 3.5,
        box.padding = 0.4,
        point.padding = 0.3,
        max.overlaps = Inf,
        min.segment.length = 0
    ) +
    scale_color_gradient(
        low = "grey85",
        high = "#F4B400"
    ) +
    labs(
        x = "SCZ MAGMA weight",
        y = "LA vs all logFC",
        color = "scDRS\ncontribution",
        title = "SCZ – LA"
    ) +
    theme_classic(base_size = 11)

pdf(file.path(plots_dir,"scdrs_SCZ_LA_gene_contributions.pdf"), width = 4.5, height = 4)
p
dev.off()

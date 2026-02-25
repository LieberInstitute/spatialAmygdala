# modified from: https://github.com/LieberInstitute/spatialdACC/blob/main/code/11_differential_expression/nnSVG_precast_DE.R
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library('EnhancedVolcano')
    library("patchwork")
})

# load spe for k=9 without WM-CC
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
spe 

# load
spe_pseudo <- readRDS(here("processed-data", "Visium", "08_marker_genes","pseudo_BS_K16_final_labels_dropped.rds"))

registration_mod <-
    registration_model(spe_pseudo, covars = NULL)

block_cor <-
    registration_block_cor(spe_pseudo, registration_model = registration_mod)

results_pairwise <-
    registration_stats_pairwise(
        spe_pseudo,
        registration_model = registration_mod,
        block_cor = block_cor,
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
    )

results_enrichment <-
    registration_stats_enrichment(
        spe_pseudo,
        block_cor = block_cor,
        covars = NULL,
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
    )

results_anova <-
    registration_stats_anova(
        spe_pseudo,
        block_cor = block_cor,
        covars = NULL,
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
    )

modeling_results <- list(
    "pairwise" = results_pairwise,
    "enrichment" = results_enrichment,
    "anova" = results_anova
)

###save modeling results list
save(
    modeling_results,
    file = here("processed-data",  "Visium", "08_marker_genes", "BS_k16_modeling_results.Rdata")
)

colData(spe)$spatialLIBD <- colData(spe)$registration_variable

sig_genes <- sig_genes_extract(
    n = 50,
    modeling_results = modeling_results,
    model_type = "enrichment",
    sce_layer = spe_pseudo
)

write.csv(sig_genes, file = here::here("processed-data",  "Visium", "08_marker_genes", "BS_k16_sig_genes_50_wBLVM.csv"), row.names = FALSE)

## For sig_genes_extract_all() to work
spe_pseudo$spatialLIBD <- spe_pseudo$layer
rownames(spe_pseudo) <- rowData(spe_pseudo)$gene_id

sig_genes_all <- sig_genes_extract_all(
    n = 50,
    modeling_results = modeling_results,
    sce_layer = spe_pseudo
)

saveRDS(sig_genes_all, file = here::here("processed-data",  "Visium", "08_marker_genes", "BS_k16_sig_genes_all_wBLVM.rds"))



# find genes that are repeated in the "gene" column
repeated_genes <- sig_genes[duplicated(sig_genes$gene),]
# top model_type test       gene     stat         pval          fdr
# 89   39 enrichment   BA      KCNG2 4.675013 1.109825e-05 3.502053e-03
# 255   5 enrichment  CoA      MOXD1 7.185029 2.524511e-10 7.966094e-07
# 270  20 enrichment  CoA    ZDHHC14 5.734806 1.497403e-07 1.232623e-04
# 271  21 enrichment  CoA AC116345.1 5.713910 1.635930e-07 1.290544e-04
# 272  22 enrichment  CoA     GALNT8 5.690561 1.805661e-07 1.367464e-04
# 276  26 enrichment  CoA     PCDH17 5.471320 4.526440e-07 2.764487e-04
# 285  35 enrichment  CoA      AIPL1 5.108634 2.000691e-06 9.238800e-04
# 287  37 enrichment  CoA      TMOD1 5.095416 2.110209e-06 9.270725e-04
# 345  45 enrichment  HPC      PRKCB 4.843528 5.753652e-06 1.078553e-03
# 359   9 enrichment  ITC  LINC00523 6.447762 6.838205e-09 1.263540e-05
# 360  10 enrichment  ITC       WFS1 6.431664 7.341117e-09 1.263540e-05
# 361  11 enrichment  ITC      MEIS2 6.296999 1.326411e-08 1.931764e-05
# 364  14 enrichment  ITC       RARB 6.102393 3.096513e-08 3.257016e-05
# 372  22 enrichment  ITC      HACD2 5.782411 1.223205e-07 7.055466e-05
# 376  26 enrichment  ITC      SYNPR 5.594454 2.705754e-07 1.280701e-04
# 389  39 enrichment  ITC      GRIK3 5.127715 1.852056e-06 5.032336e-04
# 397  47 enrichment  ITC      PDE7B 5.012129 2.947533e-06 6.599169e-04
# 399  49 enrichment  ITC    TMEM272 4.982521 3.317439e-06 7.057200e-04
# 409   9 enrichment   LA       OPTN 8.293296 1.569179e-12 3.301030e-09
# 472  22 enrichment  MeA     GABRG1 5.598466 2.660868e-07 2.238915e-04
# 594  44 enrichment Chat       LHX5 6.884376 9.793934e-10 3.582361e-07
#     gene_index     logFC    ensembl
# 89       15821 2.1250511      KCNG2
# 255       6677 2.1158483      MOXD1
# 270       6793 0.7767176    ZDHHC14
# 271       5266 4.0323914 AC116345.1
# 272      11027 1.0151088     GALNT8
# 276      12171 1.0432079     PCDH17
# 285      14538 3.1257278      AIPL1
# 287       8805 0.8771677      TMOD1
# 345      13881 1.7598219      PRKCB
# 359      12872 4.6927319  LINC00523
# 360       4336 1.5983433       WFS1
# 361      13024 2.5819874      MEIS2
# 364       3322 2.2196277       RARB
# 372       3835 1.1067957      HACD2
# 376       3646 2.3429398      SYNPR
# 389        471 1.9223588      GRIK3
# 397       6695 2.0491152      PDE7B
# 399      12156 3.0308658    TMEM272
# 409       9245 0.9506873       OPTN
# 472       4475 1.2115299     GABRG1
# 594      11800 3.4447260       LHX5


indices <- c()

indices <- append(indices, which(sig_genes$gene == "FIBCD1")) #aBA
indices <- append(indices, which(sig_genes$gene == "CNR1")) # known aBA
indices <- append(indices, which(sig_genes$gene == "PEX5L")) #known BA marker
indices <- append(indices, which(sig_genes$gene == "SLIT1")) # novel BA marker
indices <- append(indices, which(sig_genes$gene == "MOXD1")) #known BLVM.1 marker, 
indices <- append(indices, which(sig_genes$gene == "GJB3")) #novel BLVM.1 marker
indices <- append(indices, which(sig_genes$gene == "FSHB")) #known BLVM.2 marker, 
indices <- append(indices, which(sig_genes$gene == "MT3")) #novel BLVM.2 marker
indices <- append(indices, which(sig_genes$gene == "PENK")) # novel CeA marker
indices <- append(indices, which(sig_genes$gene == "GPR88")) # novel CeA marker
indices <- append(indices, which(sig_genes$gene == "PTH")) # CLA marker?
indices <- append(indices, which(sig_genes$gene == "RGS12")) # CLA marker in rodents
indices <- append(indices, which(sig_genes$gene == "PDYN")) # CoA
indices <- append(indices, which(sig_genes$gene == "GRP")) # known CoA marker
indices <- append(indices, which(sig_genes$gene == "CABP7")) # HPC
indices <- append(indices, which(sig_genes$gene == "POU3F1")) # HPC
indices <- append(indices, which(sig_genes$gene == "TSHZ1")) # known ITC marker
indices <- append(indices, which(sig_genes$gene == "FOXP2")) # known ITC marker
indices <- append(indices, which(sig_genes$gene == "SATB1")) # known LA marker
indices <- append(indices, which(sig_genes$gene == "CYP26B1")) # novel LA marker
indices <- append(indices, which(sig_genes$gene == "OTP")) # known MeA marker
indices <- append(indices, which(sig_genes$gene == "SST")) # novel MeA marker
indices <- append(indices, which(sig_genes$gene == "PNOC")) # known Chat marker
indices <- append(indices, which(sig_genes$gene == "CARTPT")) # known Chat marker
indices <- append(indices, which(sig_genes$gene == "NTS")) # known WM marker
indices <- append(indices, which(sig_genes$gene == "MBP")) # known WM marker



pdf(file = here::here("plots",  "Visium", "08_marker_genes","boxplots_BS_k16_new_wBLVM.pdf"),
    width = 8.5, height = 8)

for (i in indices) {
    layer_boxplot(
        i,
        sig_genes = sig_genes,
        short_title = TRUE,
        sce_layer = spe_pseudo,
        col_bkg_box = "grey80",
        col_bkg_point = "grey40",
        col_low_box = "violet",
        col_low_point = "darkviolet",
        col_high_box = "skyblue",
        col_high_point = "dodgerblue4",
        cex = 2,
        group_var = "layer",
        assayname = "logcounts"
    )
}

dev.off()

#volcano plots
thresh_fdr <- 0.05
thresh_logfc <- log2(1.5)
fdrs_gene_ids <- rowData(spe_pseudo)$gene_id
fdrs_gene_names <- rowData(spe_pseudo)$gene_name

df_list <- list()
plot_list <- list()

pdf(file = here::here("plots", "Visium", "08_marker_genes","volcano_BS_k16_wBLVM.pdf"),
    width = 8.5, height = 8)

for (i in unique(colData(spe_pseudo)[["layer"]])) {
    print(i)

    fdrs <- modeling_results[["enrichment"]][,paste0("fdr_", i)]
    logfc <- modeling_results[["enrichment"]][,paste0("logFC_", i)]

    # Identify significant genes (low FDR and high logFC)
    sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)

    # Number of significant genes
    print(paste("Cluster", i))
    print(table(sig))

    df_list[[i]] <- data.frame(
        gene_name = modeling_results[["enrichment"]]$gene,
        logFC = logfc,
        FDR = fdrs,
        sig = sig
    )

    p <- EnhancedVolcano(df_list[[i]],
                    lab = df_list[[i]]$gene_name,
                    pointSize = 1,
                    x = 'logFC',
                    y = 'FDR',
                    FCcutoff = 1.5,
                    pCutoff = 0.05,
                    ylab = "-log10 FDR",
                    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
                                     'FDR & Log (base 2) FC'),
                    title = paste0(i, " vs. all others"),
                    subtitle = "",
                    caption = ""
    )

    plot_list[[i]] <- p
}


dev.off()


# supp figure
png(file = here::here("plots",  "Visium", "08_marker_genes","volcano_BS_k16_wBLVM.png"),
    width = 20, height = 20, unit="in", res=300)

wrap_plots(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
           plot_list[[5]],plot_list[[6]],plot_list[[7]],
           plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],
           plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],
           nrow=4, guides="collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')

dev.off()

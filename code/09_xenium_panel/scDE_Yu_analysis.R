setwd('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("BayesSpace")
    library("RColorBrewer")
    library("ggplot2")
    library("gridExtra")
    library("patchwork")
})

# Save directories
plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", "Yu")
processed_dir = here("processed-data","09_xenium_panel")

load(here("processed-data","snRNA-seq", "yu_sce_gtf.rda"), verbose = TRUE)
sce <- sce.amy
sce
# class: SingleCellExperiment 
# dim: 31908 91699 
# metadata(0):
#     assays(2): counts logcounts
# rownames(31908): MIR1302-2HG FAM138A ... AC240274.1 AC213203.1
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(91699): AAACCCAAGTCGCCAC-1 AAACCCAGTCGAGATG-1 ...
# TTTGTTGTCACTTCTA-14 TTTGTTGTCCAAATGC-14
# colData names(14): barcode library_id ... orig_anno ident
# reducedDimNames(2): PCA UMAP
# mainExpName: RNA
# altExpNames(0):

colData(sce) <- colData(sce)[ , !duplicated(colnames(colData(sce)))]
colnames(colData(sce))
# [1] "barcode"            "library_id"         "sample"            
# [4] "n_genes_by_counts"  "total_counts"       "total_counts_mt"   
# [7] "pct_counts_mt"      "doublet_scores"     "predicted_doublets"
# [10] "celltype"           "clusters"           "space_anno"        
# [13] "orig_anno"          "ident"   

unique(sce$orig_anno)
# [1] Human_Endo NOSTRIN   Human_LAMP5 ABO      Human_Astro_1 FGFR3 
# [4] Human_PVALB ADAMTS5  Human_LAMP5 COL25A1  Human_SOX11 EBF2    
# [7] Human_Oligo_5 OPALIN Human_VGLL3 CNGB1    Human_LAMP5 COL14A1 
# [10] Human_Oligo_1 OPALIN Human_Astro_4 FGFR3  Human_Astro_3 FGFR3 
# [13] Human_CALCR LHX8     Human_HGF ESR1       Human_HGF C11orf87  
# [16] Human_HGF NPSR1      Human_DRD2 PAX6      Human_VGLL3 MEPE    
# [19] Human_PRKCD          Human_TSHZ1 CALCRL   Human_LAMP5 BDNF    
# [22] Human_VIP ABI3BP     Human_HTR3A DRD2     Human_VIP NDNF      
# [25] Human_LAMP5 NDNF     Human_Astro_2 FGFR3  Human_OPC_3 PDGFRA  
# [28] Human_DRD2 ISL1      Human_SST HGF        Human_SATB2 CALCRL  
# [31] Human_TFAP2C         Human_SATB2 IL15     Human_Oligo_2 OPALIN
# [34] Human_OPC_2 PDGFRA   Human_TSHZ1 SEMA3C   Human_Micro CTSS    
# [37] Human_SST EPYC       Human_RXFP2 RSPO2    Human_Oligo_3 OPALIN
# [40] Human_OPC_4 PDGFRA   Human_SATB2 ST8SIA2  Human_STRIP2        
# [43] Human_OPC_1 PDGFRA   Human_Oligo_4 OPALIN Human_Oligo_6 OPALIN

unique(sce$celltype)
# [1] Endothelial     ExN             Astrocyte       InN            
# [5] Oligodendrocyte OPC             Microglia     

unique(sce$space_anno)
# [1] non-neuron                BLA                      
# [3] Cortical interneuron-like PL                       
# [5] COA/MEA                   CEA                      
# [7] IA                        NLOT   


library(scran)

# find and remove MT genes
mt <- grep("^MT-", rownames(sce))
sce <- sce[-mt,]



# =========== Intercalated cell clusters ===========
# change COA/MEA to MEA in $space_anno
levels(sce$space_anno) <- c(levels(sce$space_anno), "MEA")
sce$space_anno[sce$space_anno == "COA/MEA"] <- "MEA"
unique(sce$space_anno)

sce$space_anno <- droplevels(sce$space_anno, exclude = "COA/MEA")

# Check the updated levels
levels(sce$space_anno)

# find marker genes
marker.info <- findMarkers(sce, sce$space_anno, test="binom", direction="up", lfc=2)
marker.info

# get cluster id for loop
cluster_id <- unique(levels(as.factor(sce$space_anno)))

# loop through cluster ids
for (i in 1:length(cluster_id)) {
    
    # just look at top 5 genes for one spatial cluster as a sanity check
    chosen <- marker.info[[i]]
    
    plot_dir = here("plots", "09_xenium_panel", "snRNA-seq","Yu","space_anno", cluster_id[i])
    
    #make sure folder exists, if not, make it
    if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
    }
    
    # output for keeping track with loop
    print(paste("Now printing marker genes for:", cluster_id[i]))
    
    
    # visualize mean AUC
    top10 <- c(rownames(chosen[1:10,]))
    p <- scater::plotExpression(sce, features=top10, x='space_anno', colour_by='space_anno', show_violin=TRUE, show_median=TRUE, ncol=2)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_Top10.pdf"))
    print(p)
    dev.off()
    
}

# ============ CEA vs MEA ============

sce.subset <- sce[,grepl("CEA|MEA", sce$space_anno)]
sce.subset$space_anno <- droplevels(sce.subset$space_anno, exclude = "CEA|MEA")
unique(sce.subset$space_anno)
# [1] MEA CEA
# Levels: CEA MEA

# find marker genes
marker.info <- marker.info <- findMarkers(sce.subset, sce.subset$space_anno, test="t", direction="up")
marker.info

# get cluster id for loop
cluster_id <- unique(levels(as.factor(sce.subset$space_anno)))

# loop through cluster ids
for (i in 1:length(cluster_id)) {
    
    # just look at top 5 genes for one spatial cluster as a sanity check
    chosen <- marker.info[[i]]
    
    plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", "Yu", "1v1_CEA_MEA", cluster_id[i])
    
    #make sure folder exists, if not, make it
    if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
    }
    
    # output for keeping track with loop
    print(paste("Now printing marker genes for:", cluster_id[i]))
    
    
    # visualize mean AUC
    top10 <- c(rownames(chosen[1:10,]))
    p <- scater::plotExpression(sce.subset, features=top10, x='space_anno', colour_by='space_anno', show_violin=TRUE, show_median=TRUE, ncol=2)
    
    pdf(width=7, height=7, here(plot_dir, "expressionPlots_Top10.pdf"))
    print(p)
    dev.off()
    
}




# === Volano plots ===
library(EnhancedVolcano)

marker.info <- findMarkers(sce.subset, sce.subset$space_anno, test="binom")

# Loop through each unique cluster ID
for(cluster_id in unique(sce.subset$space_anno)){
    
    # set plot dir
    plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", "Yu", "1v1_CEA_MEA", cluster_id)
    
    # Filter the all_markers data frame for the current cluster
    cluster_markers <- as.data.frame(marker.info[[cluster_id]])
    
    volcano_plot <- EnhancedVolcano(
        cluster_markers,
        lab = row.names(cluster_markers),
        x = 'summary.logFC',
        y = 'p.value',
        title = paste('Volcano plot of Cluster', cluster_id, 'markers'),
        FCcutoff = 2,
        pCutoff = 10e-200,
        drawConnectors = TRUE,
        widthConnectors = 0.75
    )
    
    
    # Optionally, save the volcano plot to a file
    ggsave(filename = here(plot_dir, paste0("volcano_plot_cluster_", cluster_id, ".png")), plot = volcano_plot)
    
}




# plot expression of SST for sce.subset
plot_dir = here("plots", "09_xenium_panel", "snRNA-seq", "Yu", "1v1_CEA_MEA", "CRH")
#make sure folder exists, if not, make it
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
}

p <- scater::plotExpression(sce.subset, features="CRH", x='space_anno', colour_by='space_anno', show_violin=TRUE, show_median=TRUE, ncol=2)

pdf(width=7, height=7, here(plot_dir, "expressionPlots_Top10.pdf"))
print(p)
dev.off()



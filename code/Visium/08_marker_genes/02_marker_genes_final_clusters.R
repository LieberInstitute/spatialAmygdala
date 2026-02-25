suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("RColorBrewer")
    library("ggplot2")
    library("patchwork")
    library("DeconvoBuddies")
})

out_dir <- here("processed-data", "Visium", "08_marker_genes", "final_clusters")

# load
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_relabel_ITC_smoothed.rds"))
spe 
# class: SpatialExperiment 
# dim: 36601 226181 
# metadata(0):
# assays(2): counts logcounts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type gene_search
# colnames(226181): AAACAAGTATCTCCCA-1_V13Y24-346_A1
#   AAACAATCTACTAGCA-1_V13Y24-346_A1 ... TTGTTTCATTAGTCTA-1_V13B23-407_C1
#   TTGTTTCCATACAACT-1_V13B23-407_C1
# colData names(59): sample_id in_tissue ... ITC_gmm2_smooth
#   BS_k16_relabel_ITC_smooth
# reducedDimNames(2): PCA PCA-HARMONY_sample
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


# relabel LA in Br6660 to CLA. $BS_k16_relabel_ITC_smooth clusters
spe$BS_k16_relabel_ITC_smooth <- as.character(spe$BS_k16_relabel_ITC_smooth)
spe$BS_k16_relabel_ITC_smooth[spe$sample_id == "Br6660" & spe$BS_k16_relabel_ITC_smooth == "LA"] <- "CLA"
spe$BS_k16_relabel_ITC_smooth <- factor(spe$BS_k16_relabel_ITC_smooth)

# rename Other to Chat
spe$BS_k16_relabel_ITC_smooth <- plyr::revalue(spe$BS_k16_relabel_ITC_smooth, c("Other" = "CHAT"))
spe$BS_k16_relabel_ITC_smooth <- plyr::revalue(spe$BS_k16_relabel_ITC_smooth, c("ITC" = "AI"))
spe$BS_k16_relabel_ITC_smooth <- plyr::revalue(spe$BS_k16_relabel_ITC_smooth, c("aBA" = "BM"))
spe$BS_k16_relabel_ITC_smooth <- plyr::revalue(spe$BS_k16_relabel_ITC_smooth, c("Meninges" = "Endothelial"))
spe$BS_k16_relabel_ITC_smooth <- plyr::revalue(spe$BS_k16_relabel_ITC_smooth, c("BLVM.1" = "PL"))
spe$BS_k16_relabel_ITC_smooth <- plyr::revalue(spe$BS_k16_relabel_ITC_smooth, c("BLVM.2" = "BL"))
spe$BS_k16_relabel_ITC_smooth <- plyr::revalue(spe$BS_k16_relabel_ITC_smooth, c("BA" = "BLD"))

# combine BLVM.1 and BLVM.2
#spe$BS_k16_relabel_ITC_smooth <- plyr::revalue(spe$BS_k16_relabel_ITC_smooth, c("BLVM.1" = "BLVM", "BLVM.2" = "BLVM"))
table(spe$BS_k16_relabel_ITC_smooth)
      #    BM         BLD          PL          BL         CeA         CLA 
      # 10363       11331       10638       49673       11678        4834 
      #   CoA         HPC          AI          LA         MeA Endothelial 
      # 12183        4774        1578       40228       13913        2611 
      #  CHAT        WM.1        WM.2 
      #  3076       29732       19569 

## Convert from character to a factor
spe$BS_k16_relabel_ITC_smooth <- as.factor(spe$BS_k16_relabel_ITC_smooth)

# Make new final name for cluster colname
spe$BS_k16_Semisupervised_wAI <- spe$BS_k16_relabel_ITC_smooth
spe$BS_k16_relabel_ITC_smooth <- NULL

colnames(colData(spe))
# 31] "exclude_overlapping"         "sum_umi_log"                
# [33] "sum_umi_outliers"            "sum_umi_z"                  
# [35] "sum_gene_log"                "sum_gene_outliers"          
# [37] "sum_gene_z"                  "expr_chrM_ratio_outliers"   
# [39] "expr_chrM_ratio_z"           "local_outliers"             
# [41] "sizeFactor"                  "BS_k10"        

# flag spots with < 50 UMIs as TRUE/FALSE
spe$low_umi <- spe$sum_umi < 75
spe$low_gene <- spe$sum_gene < 75

# drop all outliers plus low UMI
spe <- spe[, !(spe$low_umi | spe$low_gene)]
spe


# save 
saveRDS(spe, file=here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes","spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))


# ======== MARKER GENES =========

# renormalize
markers <- scran::findMarkers(
  spe,
  groups=spe$BS_k16_Semisupervised_wAI,
  test.type = c("wilcox"),
  pval.type = c("all"),
  full.stats = TRUE,
  sorted = TRUE,
  direction="up"
)

# get top 10 gene names for each cluster
top_markers <- lapply(markers, function(x) {
  top10 <- x[1:10,]
  top10_genes <- rownames(top10)
  top10_genes
})
# $`1`
#  [1] "PENK"    "SYNPR"   "TMEM272" "GAD2"    "SLC32A1" "SLC35F1" "GPR88"  
#  [8] "NLRP1"   "TSHZ1"   "GAD1"   

# save markers as csv
saveRDS(markers, file=here("processed-data","Visium","08_marker_genes", "markers_bs_final_ITC_smoothed.rds"))


# csv of top 100 markers per cluster

top100 <- lapply(markers, function(x) {
  top100 <- x[1:100,]
  top100_genes <- rownames(top100)
  top100_genes
})

# convert to data frame
top100_df <- do.call(rbind, lapply(names(top100), function(cluster) {
  data.frame(
    cluster = cluster,
    gene = top100[[cluster]],
    stringsAsFactors = FALSE
  )
}))

# write to csv  
write.csv(top100_df, file=here("processed-data","Visium","08_marker_genes", "top100_markers_bs_final_ITC_smoothed.csv"), row.names = FALSE)

# BADL = PEX5L/EDIL3/STXBP6/COL25A1
# BL = MTPN/COL25A1/NRXN2
# BLVM = LAMP5/NPTX1/STMN4/ATP2B4
# BM = NCAM2/CNR1/CCK/GABRD
# PCo = PDYN/CDH13/ESR1
# Ce = PENK/SYNPR/GAD2/TSHZ1/SST/PRKCD
# LA = CAMK2N1/TTC9B/CYP26B1/SLC30A3/ARPP19
# Me = CALB2/CALB1/ CARTPT/SLC17A6/GABRE/OTP

features <- c("PEX5L", "EDIL3", "STXBP6", "COL25A1",
                "MTPN",  "NRXN2",
                "LAMP5", "NPTX1", "STMN4", "ATP2B4",
                "NCAM2", "CNR1", "CCK", "GABRD",
                "PDYN", "CDH13", "ESR1",
                "PENK", "SYNPR", "GAD2", "SST", "PRKCD",
                "CAMK2N1", "TTC9B", "CYP26B1", "SLC30A3", "ARPP19",
                "CALB2", "CALB1", "CARTPT", "SLC17A6", "GABRE", "OTP",
                "FOXP2", "TSHZ1", "DRD1", "OPRM1", "CPNE4", "PRKG1", "SIM1", "GULP1")


# grouped heatmap
pdf(here("plots", "Visium", "08_marker_genes","final_clusters" ,"grouped_heatmap.pdf"), width = 10, height = 10)
p <- scater::plotDots(spe, features, group="BS_k16_Semisupervised_wAI", center=TRUE, scale=TRUE)
print(p)
dev.off()


# heatmap of known marker genes
known <- c("GULP1", "SATB1","COL25A1","PEX5L", "ESR1", "MOXD1","RORB","PENK","SST","FOXP2", "TSHZ1","DRD1")

# blue to red coloscale
pdf(here("plots", "Visium", "08_marker_genes","final_clusters", "amygdala_heatmap_known_markers.pdf"), width = 4, height = 3)
p <- scater::plotGroupedHeatmap(spe, known, group="BS_k16_Semisupervised_wAI", center=TRUE, scale=TRUE, treeheight_row = 0, treeheight_col = 0)
print(p)
dev.off()


novel <- c("STXBP6", "EDIL3", "NRXN2", "MTPN", "STMN4", "LAMP5", "CNR1", "NCAM2", "PENK", "NTS", "FOXP2", "TSHZ1", "CYP26B1", "TTC9B", "CALB2", "CARTPT", "PDYN", "CDH13")

# blue to red coloscale
pdf(here("plots", "Visium", "08_marker_genes","final_clusters", "amygdala_heatmap_novel_markers.pdf"), width = 4, height = 5)
scater::plotGroupedHeatmap(spe, novel, group="BS_k16_Semisupervised_wAI", center=TRUE, scale=TRUE, cluster_cols=FALSE, cluster_rows=FALSE, legend=TRUE)

dev.off()


# look for enrichment of immature neuron markers from Velmeshev et al 2019
# Immature neurons: DCX, SOX11, BCL2, EOMES
immature <- c("DCX", "SOX11", "BCL2", "EOMES")
pdf(here("plots", "Visium", "08_marker_genes","final_clusters", "amygdala_heatmap_immature_neuron_markers.pdf"), width = 4, height = 2)
scater::plotGroupedHeatmap(spe, immature, group="BS_k16_Semisupervised_wAI", cluster_cols=FALSE, cluster_rows=FALSE, legend=TRUE)
dev.off()




# ======== QC plots by cluster ========
library(scater)

plot_dir <- here("plots", "Visium", "08_marker_genes", "cluster_QC")


pal <- c( 
  AI        = "#D62728",  # strong red standout
  BM        = "#E67E22",  # burnt orange
  BLD         = "#9B59B6",  # violet
  PL     = "#f1e438ff",  # steel blue
  BL     = "#035185ff",  # sky blue
  LA         = "#F4B400",  # goldenrod
  CoA        = "#5DA5DA",  # blue-gray
  CeA        = "#197d43ff",  # emerald green
  MeA        = "#baf739ff",  # chartreuse
  HPC        = "#d6a8f8ff",  # deep royal purple (distinct from BA)
  CHAT       = "#A0522D",  # chestnut brown
  Endothelial   = "#444444",  # charcoal gray
  WM.1       = "#BBBBBB",  # light slate gray
  WM.2       = "#DDDDDD",   # mist gray
  CLA        = "#FF69B4"   # hot pink
)


# violin plots of QC metrics by final cluster
qc_metrics <- c("sum_umi", "sum_gene", "expr_chrM_ratio")

# mito ratio
png(width=15, height=5, here(plot_dir,"Violin_mito_ratio.png"), res=300, units="in")
p1 <- plotColData(spe, x="BS_k16_Semisupervised_wAI", y="expr_chrM_ratio", colour_by="BS_k16_Semisupervised_wAI") + 
    #scale_y_log10() + 
    ggtitle("Mitochondrial Percent") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none") +
    scale_color_manual(values=pal) +
    ylab("Mitochondrial Percent") 

p1
dev.off()

# sum umi
png(width=10, height=5, here(plot_dir,"Violin_sum_umi.png"), res=300, units="in")
p2 <- plotColData(spe, x="BS_k16_Semisupervised_wAI", y="sum_umi", colour_by="BS_k16_Semisupervised_wAI") + 
    scale_y_log10() + 
    ggtitle("Sum UMI") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")+
    scale_color_manual(values=pal) +
    ylab("Total UMI") 

p2
dev.off()

# sum gene
png(width=10, height=5, here(plot_dir,"Violin_sum_gene.png"), res=300, units="in")
p3 <- plotColData(spe, x="BS_k16_Semisupervised_wAI", y="sum_gene", colour_by="BS_k16_Semisupervised_wAI") + 
    scale_y_log10() + 
    ggtitle("Sum Gene") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")+
    scale_color_manual(values=pal) +
    xlab("Spatial Domains") +
    ylab("Unique Genes")

p3
dev.off()

png(width=10, height=10, here(plot_dir,"Violin_all_metrics.png"), res=300, units="in")
p1 <- p1 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank()
                 )

p2 <- p2 + theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()
                )

p1/p2/p3
dev.off()


# plot again but add A, B, C to panels. Using some other package
library(ggpubr)

png(width=10, height=10, here(plot_dir,"Violin_all_metrics_labeled.png"), res=300, units="in")
p1 <- p1 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank()
) 

p2 <- p2 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank()
)

ggarrange(p1, p2, p3,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)
dev.off()
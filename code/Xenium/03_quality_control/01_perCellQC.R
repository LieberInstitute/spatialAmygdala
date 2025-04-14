library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")

processed_dir <- here("processed-data", "Xenium", "03_quality_control")
plot_dir <- here("plots", "Xenium", "03_quality_control")

load(here("processed-data","Xenium", "02_build_spe", "spe_combined.Rdata"))
spe
# class: SpatialExperiment 
# dim: 541 1018069 
# metadata(0):
# assays(1): counts
# rownames(541): ENSG00000069431 ENSG00000151388 ...
#   DeprecatedCodeword_0344 DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(1018069): aaaadgkh-1 aaaadlfe-1 ... oihobmbk-1 oihoeehh-1
# colData names(11): cell_id transcript_counts ... sample_id brnum
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id


rowData(spe) 
# DataFrame with 541 rows and 3 columns
#                                             ID                 Symbol
#                                    <character>            <character>
# ENSG00000069431                ENSG00000069431                  ABCC9
# ENSG00000151388                ENSG00000151388               ADAMTS12
# ENSG00000145536                ENSG00000145536               ADAMTS16
# ENSG00000156140                ENSG00000156140                ADAMTS3
# ENSG00000120907                ENSG00000120907                 ADRA1A
# ...                                        ...                    ...
# DeprecatedCodeword_0196 DeprecatedCodeword_0.. DeprecatedCodeword_0..
# DeprecatedCodeword_0304 DeprecatedCodeword_0.. DeprecatedCodeword_0..
# DeprecatedCodeword_0318 DeprecatedCodeword_0.. DeprecatedCodeword_0..
# DeprecatedCodeword_0344 DeprecatedCodeword_0.. DeprecatedCodeword_0..
# DeprecatedCodeword_0373 DeprecatedCodeword_0.. DeprecatedCodeword_0..
#                                        Type
#                                 <character>
# ENSG00000069431             Gene Expression
# ENSG00000151388             Gene Expression
# ENSG00000145536             Gene Expression
# ENSG00000156140             Gene Expression
# ENSG00000120907             Gene Expression

sce <- scuttle::addPerCellQC(sce,subsets = list(Mito = which(seqnames(sce) == "chrM")))

# make scuttle subsets for rowData(sce)$Type= Gene Expression

gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe <- scuttle::addPerCellQC(spe, subsets = list(GEX = gene_expression_idx))

colnames(colData(spe))
#  [1] "cell_id"                    "transcript_counts"         
#  [3] "control_probe_counts"       "control_codeword_counts"   
#  [5] "unassigned_codeword_counts" "deprecated_codeword_counts"
#  [7] "total_counts"               "cell_area"                 
#  [9] "nucleus_area"               "sample_id"                 
# [11] "brnum"                      "sum"                       
# [13] "detected"                   "subsets_GEX_sum"           
# [15] "subsets_GEX_detected"       "subsets_GEX_percent"       
# [17] "total"      


# plot violins of GEX sum and detected
png(file = here("plots","Xenium", "03_quality_control", "GEX_sum_detected_violin.png"), width = 10, height = 5, units = "in", res = 300)
p1 <- plotColData(spe, x = "brnum", y = "subsets_GEX_sum", colour_by = "brnum") +
    ggtitle("Library size")
p2 <- plotColData(spe, x = "brnum", y = "subsets_GEX_detected", colour_by = "brnum") +
    ggtitle("Detected genes")
p1 + p2
dev.off()


# make cell and nucleus area scaling factors
spe$cell_area.sf <- spe$cell_area / median(spe$cell_area)
spe$nucleus_area.sf <- spe$nucleus_area / median(spe$nucleus_area)

# histograms of scaling factor
png(file = here("plots","Xenium", "03_quality_control", "cell_area_scaling_factors.png"), width = 10, height = 5, units = "in", res = 300)
hist(spe$cell_area.sf, breaks = 50, main = "Cell area", xlab = "Scaling factor")
dev.off()


png(file = here("plots","Xenium", "03_quality_control", "nucleus_area_scaling_factors.png"), width = 10, height = 5, units = "in", res = 300)
hist(spe$nucleus_area.sf, breaks = 50, main = "Nucleus area", xlab = "Scaling factor")
dev.off()


###### Note: Nucleus area scaling factors look most similar to what is shown in Jean Fan paper. Going to use that by default. 
# see here: https://github.com/LylaAtta123/normalization-analyses/blob/main/R/xenium.ipynb
######

# normalize the counts by the nucleus and cell area scaling factors
assay(spe, "nucleus_normcounts") <- scuttle::normalizeCounts(spe, size.factors=spe$nucleus_area.sf, transform="log", assay.type="counts")
assay(spe, "cell_normcounts") <- scuttle::normalizeCounts(spe, size.factors=spe$cell_area.sf, transform="log", assay.type="counts")

# save spe
save(spe, file = here("processed-data","Xenium", "03_quality_control", "spe_normcounts.Rdata"))


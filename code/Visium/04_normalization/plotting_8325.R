suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("spatialLIBD"))


load(here("processed-data", "Visium", "04_normalization", "spe_stitched_norm.Rdata"), verbose = TRUE)
spe

# discard out of tissue spots
spe <- spe[,colData(spe)$in_tissue]
rownames(spe) <- rowData(spe)$symbol

# exclude overlapping and NA
spe <- spe[,!colData(spe)$exclude_overlapping]
spe <- spe[,!is.na(colData(spe)$exclude_overlapping)]
spe

# Cea Marker genes
cea_markers <- c("SST", "GAD1", "GAD2","PENK")
png(here("plots", "Visium", "04_normalization", "cea_markers.png"), width=10, height=10, units="in", res=300)
vis_gene(spe, 
    geneid = cea_markers, 
    sampleid = "Br8325",
    is_stitched=TRUE,
    multi_gene_method="pca"
    )
dev.off()   


png(here("plots", "Visium", "04_normalization", "cea_markers_penk.png"), width=7, height=7, units="in", res=300)
vis_gene(spe, 
    geneid = "PENK", 
    sampleid = "Br8325",
    is_stitched=TRUE,
    multi_gene_method="pca"
    )
dev.off()   


# Mea Marker genes
mea_markers <- c("SLC17A6","SIM1","CARTPT")
png(here("plots", "Visium", "04_normalization", "mea_markers.png"), width=10, height=10, units="in", res=300)
vis_gene(spe, 
    geneid = mea_markers, 
    sampleid = "Br8325",
    is_stitched=TRUE,
    multi_gene_method="pca"
    )
dev.off()   

png(here("plots", "Visium", "04_normalization", "mea_markers_slc17a6.png"), width=7, height=7, units="in", res=300)
vis_gene(spe, 
    geneid = "SLC17A6", 
    sampleid = "Br8325",
    is_stitched=TRUE,
    multi_gene_method="pca"
    )
dev.off()   


# ITC Marker genes
itc_markers <- c("TSHZ1","MEIS2","DRD1", "CPNE4","SHISA9")
png(here("plots", "Visium", "04_normalization", "itc_markers.png"), width=10, height=10, units="in", res=300)
vis_gene(spe, 
    geneid = itc_markers, 
    sampleid = "Br8325",
    is_stitched=TRUE,
    multi_gene_method="pca"
    )
dev.off()   


png(here("plots", "Visium", "04_normalization", "itc_markers_tshzs1.png"), width=10, height=10, units="in", res=300)
vis_gene(spe, 
    geneid = "TSHZ1", 
    sampleid = "Br8325",
    is_stitched=TRUE,
    multi_gene_method="pca"
    )
dev.off()   



png(here("plots", "Visium", "04_normalization", "SLC17A7_markers.png"), width=7, height=7, units="in", res=300)
vis_gene(spe, 
    geneid = "SLC17A7", 
    sampleid = "Br8325",
    is_stitched=TRUE,
    multi_gene_method="pca"
    )
dev.off()   


png(here("plots", "Visium", "04_normalization", "CARTPT_markers.png"), width=7, height=7, units="in", res=300)
vis_gene(spe, 
    geneid = "CARTPT", 
    sampleid = "Br8325",
    is_stitched=TRUE,
    multi_gene_method="pca"
    )
dev.off()   


png(here("plots", "Visium", "04_normalization", "SIM1_markers.png"), width=7, height=7, units="in", res=300)
vis_gene(spe, 
    geneid = "SIM1", 
    sampleid = "Br8325",
    is_stitched=TRUE,
    multi_gene_method="pca"
    )
dev.off()   







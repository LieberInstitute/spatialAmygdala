library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(scater)
library(patchwork)
library(ggpubr)


# install Harriet's forked repo with edge detection methods
devtools::load_all("~/forks/SpotSweeper-fork")
library(SpotSweeper)


# Save directories
plot_dir = here("plots", "Visium", "03_qc_metrics", "stitched")
processed_dir = here("processed-data","Visium","02_build_spe")

# save combined object as rds
load(here(processed_dir, "spe-1st.Rdata"))
spe 

# ========= Add QC metrics to spe ==========
colnames(colData(spe))

# QC metrics are already there, so we won't add them again


# =============== Calculate QC Metrics =================
#spe <- spe[,spe$in_tissue]

# spe$slide <- gsub("_.*", "", spe$capture_area)
# unique(spe$slide)

# ======= SpotSweeper ======
# library size
spe <- localOutliers(spe, metric="sum_umi",direction="lower", log=TRUE)

# unique genes
spe <- localOutliers(spe, metric="sum_gene", direction="lower", log=TRUE)

# mitochondrial percent
spe <- localOutliers(spe, metric="expr_chrM_ratio", direction="higher", log=FALSE)



# edge dry spots
test <- detectEdgeDryspots(spe, 
                        qc_metric = "sum_gene",
                        samples = "sample_id", 
                        mad_threshold = 3,
                        edge_threshold = 0.75,
                        shifted = FALSE,
                        #batch_var = "both",
                        name = "edge_dryspot")


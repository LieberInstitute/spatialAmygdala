library("here")
library("SingleCellExperiment")
library("scran")
library("scater")

processed_data <- here("processed-data", "Visium", "06_batch_correction")
plot_dir <- here("plots", "Visium", "06_batch_correction")

# load data with PCA embeddings
load(here(processed_data, "spe_harmony.Rdata"))
spe


# ======== plot PCA before and after Harmony correction ========

# plot PCA before Harmony correction
png(here(plot_dir, "PCA_before_harmony.png"), width = 2000, height = 1700, res = 300)
plotReducedDim(spe, ncomponents=4, dimred="PCA", scattermore=TRUE, colour_by="sample_id")
dev.off()

# plot PCA after Harmony correction
png(here(plot_dir, "PCA_after_harmony_sample-slide.png"), width = 2000, height = 1700, res = 300)
plotReducedDim(spe, ncomponents=4, dimred="PCA-HARMONY_sample_slide", scattermore=TRUE, colour_by="sample_id")
dev.off()

# plot PCA after Harmony correction (sample only)
png(here(plot_dir, "PCA_after_harmony_sample.png"), width = 2000, height = 1700, res = 300)
plotReducedDim(spe, ncomponents=4, dimred="PCA-HARMONY_sample", scattermore=TRUE, colour_by="sample_id")
dev.off()



# ====== Get explanatory PCs ======

#  uncorrected
expPCs <- getExplanatoryPCs(spe, dimred= "PCA",  
                                variables = c(
                                        "sum_umi",
                                        "expr_chrM_ratio",
                                        "sample_id",
                                        "capture_area",
                                        "slide_id"
                                    )
    )

png(here(plot_dir, "explanatory_PCs_PCA.png"), width = 2000, height = 1700, res = 300)
plotExplanatoryPCs(expPCs)
dev.off()


#  correct - sample + slide
expPCs <- getExplanatoryPCs(spe, dimred= "PCA-HARMONY_sample_slide",  
                            variables = c(
                                    "sum_umi",
                                    "expr_chrM_ratio",
                                    "sample_id",
                                    "capture_area",
                                    "slide_id"
                                )
    )

png(here(plot_dir, "explanatory_PCs_PCA-HARMONY_sample_slide.png"), width = 2000, height = 1700, res = 300)
plotExplanatoryPCs(expPCs)
dev.off()


#  correct - sample only
expPCs <- getExplanatoryPCs(spe, dimred= "PCA-HARMONY_sample",  
                            variables = c(
                                    "sum_umi",
                                    "expr_chrM_ratio",
                                    "sample_id",
                                    "capture_area",
                                    "slide_id"
                                )
    )

png(here(plot_dir, "explanatory_PCs_PCA-HARMONY_sample.png"), width = 2000, height = 1700, res = 300)
plotExplanatoryPCs(expPCs)
dev.off()


# ====== Explanatory variables =====

png(here(plot_dir, "explanatory_Vars.png"), width = 2000, height = 1700, res = 300)
plotExplanatoryVariables(spe, variables = c(
                                                    "sum_umi",
                                                    "expr_chrM_ratio",
                                                    "sample_id",
                                                    "capture_area",
                                                    "slide_id"
                                                )
                                    )
dev.off()
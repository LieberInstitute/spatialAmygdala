suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scry"))
suppressPackageStartupMessages(library("BiocSingular"))
suppressPackageStartupMessages(library("PCAtools"))
suppressPackageStartupMessages(library("patchwork"))

load(here("processed-data","04_normalization","spe_stitched_norm.Rdata"))

# normalize data


set.seed(195)
message("running PCA - ", Sys.time())

# get only top 1000 HVGs
dec <- modelGeneVar(spe)
chosen <- getTopHVGs(dec, n=4000)

spe <- scater::runPCA(spe, 
                      subset_row=chosen,
                      ncomponents = 50,
                      exprs_values='logcounts',
                      scale = TRUE, name = "PCA",
                      BSPARAM = BiocSingular::RandomParam())

#pca plots
png(here("plots", "06_dim_reduction", "PCA_sample_id_stiched.png"), width=5, height=5, units="in", res=300)
plotReducedDim(spe, dimred = "PCA", colour_by = "sample_id")
dev.off()

# make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow

png(here("plots", "06_dim_reduction", "pca_var_explained_stitched.png"), width=5, height=5, units="in", res=300)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()


# save
save(spe, file = here::here("processed-data", "06_dim_reduction", "spe_stitched_pca.Rdata"))


suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scry"))
suppressPackageStartupMessages(library("BiocSingular"))
suppressPackageStartupMessages(library("PCAtools"))
suppressPackageStartupMessages(library("patchwork"))

load(here("processed-data","03_qc_metrics","spe_stitched_local_outliers.Rdata"))

#option 1: poisson deviance feature selection, then GLM PCA
set.seed(8)
spe$sample_id <- factor(spe$sample_id)
spe <- devianceFeatureSelection(spe, assay = "counts", fam = "poisson", sorted = T)
spe

png(here("plots", "06_dim_reduction", "poisson_deviance.pdf"), width=5, height=5, units="in", res=300)
plot(sort(rowData(spe)$poisson_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylim=c(0,1000000),
     ylab = "poisson deviance", main = "Feature Selection with Deviance"
)
abline(v = 1000, lty = 2, col = "purple")
abline(v = 2000, lty = 2, col = "red")
abline(v = 2500, lty = 2, col = "pink")
abline(v = 3000, lty = 2, col = "blue")
abline(v = 4000, lty = 2, col = "green")
abline(v = 5000, lty = 2, col = "black")
dev.off()

#glm pca
hdg <- rownames(counts(spe))[1:1000]

set.seed(9)
message("running nullResiduals - ", Sys.time())
res <- spe[rownames(counts(spe)) %in% hdg,]
res <- nullResiduals(res,
                     fam = "poisson",
                     type = "pearson",
                     assay = "counts",
                      batch = res$sample_id
)

set.seed(10)
message("running PCA - ", Sys.time())

res <- scater::runPCA(res, ncomponents = 50,
                      exprs_values='poisson_pearson_residuals',
                      scale = TRUE, name = "pp-GLM-PCA")

#pca plots
png(here("plots", "06_dim_reduction", "sel_poisson_pearson_GLM_PCA_brnum.pdf"), width=5, height=5, units="in", res=300)
plotReducedDim(res, dimred = "pp-GLM-PCA", colour_by = "brnum")
dev.off()

png(here("plots", "06_dim_reduction", "sel_poisson_pearson_GLM_PCA_sample_id.pdf"), width=5, height=5, units="in", res=300)
plotReducedDim(res, dimred = "pp-GLM-PCA", colour_by = "sample_id")
dev.off()

# make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(res, "pp-GLM-PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow

png(here("plots", "06_dim_reduction", "sel_poisson_pearson_var_explained.pdf"), width=5, height=5, units="in", res=300)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

# === add to Reduced Dims and save ===
reducedDim(spe,'pp-GLM-PCA') <- reducedDim(res,'pp-GLM-PCA')


#option 2: binomial deviance feature selection, then GLM PCA

set.seed(8)
spe <- devianceFeatureSelection(spe, assay = "counts", fam = "binomial", sorted = T) 

pdf(here("plots", "06_preprocessing", "binomial_deviance.pdf"))
plot(sort(rowData(spe)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylim=c(0,1000000),
     ylab = "binomial deviance", main = "Feature Selection with Deviance"
)
abline(v = 1000, lty = 2, col = "purple")
abline(v = 2000, lty = 2, col = "red")
abline(v = 2500, lty = 2, col = "pink")
abline(v = 3000, lty = 2, col = "blue")
abline(v = 4000, lty = 2, col = "green")
abline(v = 5000, lty = 2, col = "black")
dev.off()

hdg <- rownames(counts(spe))[1:1000]

set.seed(9)
message("running nullResiduals - ", Sys.time())
res <- spe[rownames(counts(spe)) %in% hdg,]
res <- nullResiduals(res,
                     fam = "binomial",
                     type = "deviance",
                     assay='counts'
)

set.seed(10)
res <- scater::runPCA(res, ncomponents = 50,
                      exprs_values='binomial_deviance_residuals',
                      scale = TRUE, name = "bd-GLM-PCA")

percent.var <- attr(reducedDim(res, "bd-GLM-PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow

plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")

pdf(here("plots", "06_preprocessing", "binomial_deviance_GLM_PCA_brnum.pdf"))
plotReducedDim(res, dimred = "bd-GLM-PCA", colour_by = "brnum")
dev.off()

pdf(here("plots", "06_preprocessing", "binomial_deviance_GLM_PCA_sample_id.pdf"))
plotReducedDim(res, dimred = "bd-GLM-PCA", colour_by = "sample_id")
dev.off()

# === add to Reduced Dims and save ===
reducedDim(spe,'bd-GLM-PCA') <- reducedDim(res,'bd-GLM-PCA')

save(spe, file = here::here("processed-data", "06_dim_reduction", "spe_GLM-PCA.Rdata"))
library(SpatialExperiment)
library(RcppML)
library(here)
library(Matrix)


plot_dir <- here("plots","10_NMF", "NMF_LIBD")
processed_dir <- here("processed-data", "snRNAseq")

# get toy brain data
sce <- readRDS(here(processed_dir, "sce.human_all_genes.rds"))
sce

# get logcounts
logcounts <- logcounts(sce)


# run NMF
print("Starting NMF!")
start_time <- Sys.time()
x<-nmf(logcounts,
       100,
       tol = 1e-06,
       maxit = 1000,
       verbose = TRUE,
       seed = 1512,
       L1 = c(0, 0),
       mask_zeros = FALSE,
       diag = TRUE,
       nonneg = TRUE
)
end_time <- Sys.time()
print(end_time - start_time)
# Time difference of 19.1871 mins

save(x,file=here("processed-data", "Visium","10_NMF", "NMF_LIBD","RcppML_NMF_LIBD_human_allGenes.rda"))
library(SpatialExperiment)
library(scater)
library(RcppML)
library(ggspavis)
library(here)
library(scRNAseq)
library(Matrix)
library(scran)
library(scuttle)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(bluster)
library(patchwork)
library(cowplot)


plot_dir <- here("plots","10_NMF", "NMF_Yu")
processed_dir <- here("processed-data", "snRNA-seq")

# get toy brain data
load(here(processed_dir, "yu_sce_gtf.rda"))
sce.amy

# get logcounts
logcounts <- logcounts(sce.amy)


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

save(x,file=here("processed-data","10_NMF", "NMF_Yu","RcppML_NMF_Yu.rda"))
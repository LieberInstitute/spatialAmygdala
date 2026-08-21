#!/usr/bin/env Rscript
## 07b: edgeR paired pseudobulk, neuronal vs non-neuronal -> shrunken logFC
suppressPackageStartupMessages({library(edgeR)})
OUT <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond"

cnt <- read.delim(file.path(OUT,"pseudobulk_counts.tsv"), row.names=1, check.names=FALSE)
sn  <- colnames(cnt)
donor <- factor(sub("__.*","",sn)); comp <- factor(sub(".*__","",sn), levels=c("nonneuronal","neuronal"))
cat("samples:", ncol(cnt), " donors:", nlevels(donor), " design: ~donor + compartment\n")

y <- DGEList(counts=as.matrix(cnt), group=comp)
keep <- filterByExpr(y, group=comp)                    # filter on the biological contrast
cat("genes:", nrow(y), " kept by filterByExpr:", sum(keep), "\n")
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")

design <- model.matrix(~ donor + comp)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
coefname <- colnames(design)[ncol(design)]
cat("testing coefficient:", coefname, "(positive = neuronal-biased)\n")

## shrunken logFC via predFC -- stabilises low-count genes
lfc_shrunk <- predFC(y, design, prior.count=5)[, ncol(design)]
res <- glmQLFTest(fit, coef=ncol(design))
tt  <- topTags(res, n=Inf, sort.by="none")$table

out <- data.frame(gene=rownames(y), delta=as.numeric(lfc_shrunk),
                  logFC_raw=tt$logFC, logCPM=tt$logCPM, FDR=tt$FDR,
                  stringsAsFactors=FALSE)
write.table(out, file.path(OUT,"gene_neuronal_axis.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cat("\n--- delta quantiles ---\n"); print(round(quantile(out$delta, c(0,.01,.25,.5,.75,.99,1)),3))

WM <- c("MBP","PLP1","MOBP","MAG","MOG","OLIG1","OLIG2","CLDN11","CNP","SOX10")
GM <- c("RBFOX3","SNAP25","SYT1","SLC17A7","GAD1","GAD2","GRIN1","SYN1","MEG3")
EN <- c("CLDN5","FLT1","VWF","PECAM1","COL4A1")
show <- function(lab, g){ s <- out[out$gene %in% g, c("gene","delta","logCPM","FDR")]
  cat("\n---", lab, "(delta<0 = non-neuronal) ---\n"); print(s[order(s$delta),], row.names=FALSE) }
show("myelin / oligodendrocyte", WM); show("neuronal", GM); show("endothelial", EN)

## donor consistency: per-donor logCPM difference
lcpm <- cpm(y, log=TRUE, prior.count=2)
d <- sapply(levels(donor), function(dd) {
  i <- which(donor==dd & comp=="neuronal"); j <- which(donor==dd & comp=="nonneuronal")
  if (length(i)&&length(j)) lcpm[,i]-lcpm[,j] else rep(NA_real_, nrow(lcpm)) })
cat("\n--- donor consistency: correlation of per-donor neuronal-vs-non contrasts ---\n")
print(round(cor(d, use="pairwise.complete.obs"), 2))
cat("\nmean off-diagonal r:", round(mean(cor(d,use='pairwise.complete.obs')[upper.tri(diag(ncol(d)))]),3), "\n")

cat("\n--- extreme deltas: are they low-abundance? ---\n")
ex <- out[order(-abs(out$delta)),][1:15, c("gene","delta","logCPM")]
print(ex, row.names=FALSE)
cat("median logCPM all genes:", round(median(out$logCPM),2),
    " | median logCPM top-15 |delta|:", round(median(ex$logCPM),2), "\n")
cat("\nDONE\n")

library(DESeq2)
set.seed(1)
dds <- makeExampleDESeqDataSet(n=1000, m=8, betaSD=rep(0:1,c(900,100)), interceptMean=8,
#                               dispMeanRel=function(x) .1/x + 10^(rnorm(1000,-4,2)))
                               dispMeanRel=function(x) .1/x + 10^(rnorm(1000,-4,.5)))

keep <- rowSums(counts(dds) >= 10) >= 4
table(keep)
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)

dds <- DESeq(dds)
res <- results(dds)
lfc <- lfcShrink(dds, coef=2, quiet=TRUE)

counts <- t(counts(dds))
coldata <- as.data.frame(colData(dds)[,1,drop=FALSE])

write.csv(counts, file="counts_df.csv", quote=FALSE)
write.csv(coldata, file="clinical_df.csv", quote=FALSE)

system("python3.9 test_pydeseq2.py")

scan("fitted_disp_prior_var.txt")
attr(dispersionFunction(dds), "dispPriorVar")

read.csv("fitted_trend_coefs.csv", header=TRUE, row.names=1)
attr(dispersionFunction(dds), "coefficients")

disps <- read.csv("fitted_disp_table.csv", row.names=1)
mcols(dds)[,c("dispGeneEst","dispFit","dispMAP")]
head(disps)

plotDispEsts(dds)
plot(res$baseMean, disps$gene, log="xy", cex=.2)
points(res$baseMean, disps$fitted, col="red")
points(res$baseMean, disps$map, col="dodgerblue")

pyres <- read.csv("fitted_results_table.csv", row.names=1)
res

pylfc <- read.csv("fitted_shrunken_lfc.csv", row.names=1)
lfc

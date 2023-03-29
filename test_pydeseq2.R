library(DESeq2)
set.seed(1)
dds <- makeExampleDESeqDataSet(n=1000, m=12, betaSD=rep(0:1,c(400,100)),
                               interceptMean=8, interceptSD=1,
                               dispMeanRel=function(x) 5/x + 10^(rnorm(1000,-2,1)))
                               
keep <- rowSums(counts(dds) >= 10) >= 6
table(keep)
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds, alphaInit=1)
dds <- estimateDispersionsFit(dds)
dds <- estimateDispersionsMAP(dds)
plotDispEsts(dds)

dds <- nbinomWaldTest(dds)
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
points(res$baseMean, disps$fitted, col="red", cex=.2)
points(res$baseMean, disps$map, col="dodgerblue", cex=.2)

plot(mcols(dds)$dispGeneEst, disps$gene, log="xy");abline(0,1)
plot(mcols(dds)$dispMAP, disps$map, log="xy");abline(0,1)

pyres <- read.csv("fitted_results_table.csv", row.names=1)
res
head(pyres)
plot(res$stat, pyres$stat); abline(0,1)

pylfc <- read.csv("fitted_shrunken_lfc.csv", row.names=1)
lfc
head(pylfc)
plot(lfc$log2FoldChange, pylfc$log2FoldChange); abline(0,1)

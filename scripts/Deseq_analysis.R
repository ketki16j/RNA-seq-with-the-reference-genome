# Load DESeq2 library
library("DESeq2")
# Set the working directory
directory <- "~/Documents/counts2"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "Gmax_DESeq2"

sampleFiles<- c("SRR391535-output_basename.counts",
                "SRR391536-output_basename.counts",
                "SRR391537-output_basename.counts",
                "SRR391538-output_basename.counts",
                "SRR391539-output_basename.counts",
                "SRR391541-output_basename.counts")
sampleNames <- c("Leaf tissue 1","Leaf tissue 2","Leaf tissue 3","Leaf tissue 4","Leaf tissue 5","Leaf tissue 6")
sampleCondition <- c("ambient","ambient","elevated","elevated","elevated","ambient")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

treatments = c("ambient","elevated")

library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

#guts
dds <- DESeq(ddsHTSeq)
res <- results(dds)

# copied from: https://benchtobioinformatics.wordpress.com/category/dexseq/
# order results by padj value (most significant to least)
res= subset(res, padj<0.05)
res <- res[order(res$padj),]
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))

# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"))

# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.05,
             cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)
write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean <- results(ddsClean)
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

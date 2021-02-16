# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.12")
# BiocManager::install("DESeq2")

library("DESeq2")
#1 read data
setwd("~/Dropbox/course/BIEN3320/Tutorial/T03/")
countMatrix <- as.matrix(read.csv("count.csv",sep=",",row.names="gene"))
info <- read.csv("info.csv", row.names=1)
info$condition <- factor(info$condition)

# check the sample order
all(rownames(info) == colnames(countMatrix))

###2 DESeq2
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = info,
                              design = ~ condition)
dds

#3 pre-filtering
dds <- dds[rowSums(counts(dds)) > 0,]

#4 factor levels
dds$condition <- factor(dds$condition, levels = c("control","tumor"))

#5 Differential expression analysis
dds <- DESeq(dds)
results(dds)
comparison=resultsNames(dds)[2]
res <- results(dds, name=comparison)  #condition_tumor_vs_control; log2(tumor/control).  
print(comparison)
##6 save to csv file
res<-as.data.frame(res[complete.cases(res),])
write.csv(res,'./deseq2.csv')

counts(dds, normalized=TRUE)
sizeFactors(dds)

#Volcano Plot to show the differential expressed genes
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# Add colored points: blue if padj<0.05 & log2FoldChange<0, red if pvalue<.05 & log2FoldChange>0)
with(subset(res, pvalue<.05 & log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", cex = 2))
with(subset(res, pvalue<.05 & log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex = 2))


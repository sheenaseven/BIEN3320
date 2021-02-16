# Tutorial 03: featureCounts and DESeq2


## Ownership
[Wang Lab at HKUST](http://wang-lab.ust.hk)


## Data Illustration
The raw sequencing data of these patients will be only used for tutorials. All data will be deleted after this semester’s BIEN3320/LIFS4320 course.


## Introduction
• **featureCounts** is a read summarization program suitable for counting reads generated from either RNA or DNA sequencing experiments. [Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930.](https://doi.org/10.1093/bioinformatics/btt656)

• **DESeq2** is a method for differential analysis of count data. [Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), 1-21.](https://doi.org/10.1186/s13059-014-0550-8) The DESeq2 package is available at (http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html).

## Prerequisites of **DESeq2**

• A count matrix (DESeq2/Input/count.csv)

• A sample information table (DESeq2/Input/info.csv)


### step 1 

Open your RStudio and create a R script.


### step 2 

Find the directory of featureCounts (eg. yourpath/subread-2.0.1-MacOS-x86_64/bin/featureCounts)

### step 3

Run featureCounts

```
yourpath/subread-2.0.1-MacOS-x86_64/bin/featureCounts -p -T 12 -B -t exon -g gene_name \
-a yourpath/featureCounts/Input/gtf/chr21.hg37.gtf \
-o outputpath/Control1.txt yourpath/featureCounts/Input/bam/Control1.bam
```

`Note`: you need to change the **path** to your path. 

## Prerequisites of **DESeq2**

• A count matrix (DESeq2/Input/count.csv)

• A sample information table (DESeq2/Input/info.csv)


### step 1 

Install RStudio, open your RStudio and create a R script (File->New File->R script).

### step 2 

Install the “DESeq2” package. 

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library("DESeq2")
```

### step 3

Prepare your data: A count matrix and A sample information table

```
setwd("~/Dropbox/course/BIEN3320/Tutorial/T03/")

countMatrix <- as.matrix(read.csv("count.csv",sep=",",row.names="gene"))
info <- read.csv("info.csv", row.names=1)
info$condition <- factor(info$condition)

#The columns of the count matrix and the rows of the information table are in the same order.
all(rownames(info) == colnames(countMatrix)) 
```

### step 4

Run DESeq2: two-group comparison

```
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = info,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 0,]
dds$condition <- factor(dds$condition, levels = c("control","tumor"))
dds <- DESeq(dds)
results(dds)
comparison=resultsNames(dds)[2]
res <- results(dds, name=comparison)
print(comparison)
res<-as.data.frame(res)
write.csv(res,'./deseq2.csv')
```

### step 5

Volcano Plot to show the differentially expressed genes

```
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# Add colored points: blue if padj<0.05 & log2FoldChange<0, red if pvalue<.05 & log2FoldChange>0)
with(subset(res, pvalue<.05 & log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, pvalue<.05 & log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
```

<div align=center><img width="780" height="548" src="https://github.com/sheenaseven/BIEN3320/blob/main/DESeq2/Output/Volcano.plot.png"/></div>

15 Feb 2021
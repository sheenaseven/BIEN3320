# Tutorial 03: featureCounts and DESeq2


## Ownership
[Wang Lab at HKUST](http://wang-lab.ust.hk)


## Data Illustration
The raw sequencing data of these patients will be only used for tutorials. All data will be deleted after this semester’s BIEN3320/LIFS4320 course.


## Introduction
• **featureCounts** is a read summarization program suitable for counting reads generated from either RNA or DNA sequencing experiments. [Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930.](https://doi.org/10.1093/bioinformatics/btt656)

• **DESeq2** is a method for differential analysis of count data. [Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), 1-21.](https://doi.org/10.1186/s13059-014-0550-8) The DESeq2 package is available at (http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html).

## Prerequisites of **featureCounts**

• A bam file (featureCounts/Input/bam/Control1.bam)

• A GTF file (featureCounts/Input/gtf/chr21.hg37.gtf)


### step 1 

Download the package (https://sourceforge.net/projects/subread/files/subread-2.0.1/) and unzip the package.

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




15 Feb 2021
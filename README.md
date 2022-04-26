# m6AexpressBHM
m6Aexpress-BHM: Predicting m6A regulation of gene expression in multiple-groups context by a Bayesian Hierarchical Mixture model
# Installation Instructions
The m6AexpressBHM package is supported by R 4.1.0 or newer versions. First, you need to install the exomePeak2 package for m6A peak calling:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("exomePeak2")
```
To obtain the longest transcript, user shoul install *m6ALogisticModel* R package
```r
devtools::install_github("ZhenWei10/m6ALogisticModel")
```
# Usage Example
## For each concerned group (e.g., each viral infection; one time-point bacterial infected), we first did the peak calling by *exomePeak2* R package.
```r
library(exomePeak2)
library(m6ALogisticModel)
f1 <- "./group1_IP1.bam"
f2 <- "./group1_IP2.bam"
f3 <- "./group1_IP3.bam"
f4 <- "./group1_Input1.bam"
f5 <- "./group1_Input2.bam"
f6 <- "./group1_Input3.bam"
group1_IP_BAM <- c(f1,f2,f3)
group1_INPUT_BAM <- c(f4,f5,f6)
###peak calling for group1

fa <- "./group2_IP1.bam"
fb <- "./group2_IP2.bam"
fc <- "./group2_IP3.bam"
fd <- "./group2_Input1.bam"
fe <- "./group2_Input2.bam"
ff <- "./group2_Input3.bam"
group2_IP_BAM <-  c(fa,fb,fc)
group2_INPUT_BAM <- c(fd,fe,ff)
##peak calling for group2
```





```

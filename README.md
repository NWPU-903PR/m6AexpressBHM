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

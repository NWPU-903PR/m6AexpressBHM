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
GENE_ANNO_GTF = "./hg19_GTF/genes.gtf"

f1 <- "./group1_IP1.bam"
f2 <- "./group1_IP2.bam"
f3 <- "./group1_IP3.bam"
f4 <- "./group1_Input1.bam"
f5 <- "./group1_Input2.bam"
f6 <- "./group1_Input3.bam"
group1_IP_BAM <- c(f1,f2,f3)
group1_INPUT_BAM <- c(f4,f5,f6)
###peak calling for group1
group1_peak_calling <- peak_calling(IP_BAM=group1_IP_BAM,
                                    INPUT_BAM=group1_INPUT_BAM,
                                    GENE_ANNO_GTF=GENE_ANNO_GTF,
                                    paired_end = FALSE,
                                    Genome = "hg19"
                                    output_dir="./Group1_peakcalling")
                                                           
###Group2 MeRIP-seq data
fa <- "./group2_IP1.bam"
fb <- "./group2_IP2.bam"
fc <- "./group2_IP3.bam"
fd <- "./group2_Input1.bam"
fe <- "./group2_Input2.bam"
ff <- "./group2_Input3.bam"
group2_IP_BAM <-  c(fa,fb,fc)
group2_INPUT_BAM <- c(fd,fe,ff)
##peak calling for group2
group2_peak_calling <- peak_calling(IP_BAM=group2_IP_BAM,
                                    INPUT_BAM=group2_INPUT_BAM,
                                    GENE_ANNO_GTF=GENE_ANNO_GTF,
                                    paired_end = FALSE,
                                    Genome = "hg19"
                                    output_dir="./Group2_peakcalling")
                                
###Group2 MeRIP-seq data
fA <- "./group3_IP1.bam"
fB <- "./group3_IP2.bam"
fC <- "./group3_IP3.bam"
fD <- "./group3_Input1.bam"
fE <- "./group3_Input2.bam"
fF <- "./group3_Input3.bam"
group2_IP_BAM <-  c(fA,fB,fC)
group2_INPUT_BAM <- c(fD,fE,fF)
##peak calling for group3
group3_peak_calling <- peak_calling(IP_BAM=group3_IP_BAM,
                                    INPUT_BAM=group3_INPUT_BAM,
                                    GENE_ANNO_GTF=GENE_ANNO_GTF,
                                    paired_end = FALSE,
                                    Genome = "hg19"
                                    output_dir="./Group3_peakcalling")
                                    
##Combine peak sites information and reads count of each peak together 
Group1_peakinfor <- obtain_peakinfor(peak_infor_dir="./Group1_peakcalling")
Group2_peakinfor <- obtain_peakinfor(peak_infor_dir="./Group2_peakcalling")
Group3_peakinfor <- obtain_peakinfor(peak_infor_dir="./Group3_peakcalling")  
## To avoid isoform ambiguity, we selected m6A peak sites mapped to the longest transcript of each gene
Group1_mappeak_LTX <- map_peak_longTX(filepath="./Group1_peakcalling"
                                     annotation_file=GENE_ANNO_GTF,
                                     peak_sites_infor=Group1_peakinfor)
                                     
Group2_mappeak_LTX <- map_peak_longTX(filepath="./Group2_peakcalling"
                                     annotation_file=GENE_ANNO_GTF,
                                     peak_sites_infor=Group2_peakinfor)

Group3_mappeak_LTX <- map_peak_longTX(filepath="./Group3_peakcalling"
                                     annotation_file=GENE_ANNO_GTF,
                                     peak_sites_infor=Group3_peakinfor)
```
## Quantify methylation for each gene in each concerned group
```r
## Find the peak center,which was used to quantify the distance to stop codon
Group1_peakcenter <- findpeakcenter(annotation_file=GENE_ANNO_GTF,maplongTX_peak=Group1_mappeak_LTX)
Group2_peakcenter <- findpeakcenter(annotation_file=GENE_ANNO_GTF,maplongTX_peak=Group2_mappeak_LTX)
Group3_peakcenter <- findpeakcenter(annotation_file=GENE_ANNO_GTF,maplongTX_peak=Group3_mappeak_LTX)
## Obtain the distance of peak center to stop condon
Group1_peakcenter2stopcondon <- dist_stopcodon(target_peakcenter=Group1_peakcenter,annotation_file=GENE_ANNO_GTF)
Group2_peakcenter2stopcondon <- dist_stopcodon(target_peakcenter=Group2_peakcenter,annotation_file=GENE_ANNO_GTF)
Group3_peakcenter2stopcondon <- dist_stopcodon(target_peakcenter=Group3_peakcenter,annotation_file=GENE_ANNO_GTF)
## Quantify methylation level weighted by the distance to stop condon
### load library size factor for each group sample
load("./Group1_peakcalling/size_factor.Rdata")
Group1_sizefactor <- size_factor
load("./Group2_peakcalling/size_factor.Rdata")
Group2_sizefactor <- size_factor
load("./Group3_peakcalling/size_factor.Rdata")
Group3_sizefactor <- size_factor
Group1_methylevel_distdecay  <- gene_methy_level_distdecay()
```






```

##MeRIP-seq data alignment
MeRIP_Seq_Alignment <- scanMeripBAM(
                         bam_ip = IP_BAM,
                         bam_input = INPUT_BAM,
                         paired_end = FALSE)

##Perform peak calling on exons
###If either genome or bsgenome arguments are provided, the GC content bias correction will be performed while peak calling.
SummarizedExomePeaks <- exomePeakCalling(merip_bams = MeRIP_Seq_Alignment,
                                         gff_dir = GENE_ANNO_GTF,
                                         genome = "hg19") 

##Estimate correction factors
###the normalization methods used will be controlled by the functions estimateSeqDepth() and normalizeGC(). 
####User can combine multiple normalization approaches, such as adding the conditional quantile normalization introduced in
SummarizedExomePeaks <- estimateSeqDepth(SummarizedExomePeaks)
SummarizedExomePeaks <- normalizeGC(SummarizedExomePeaks)
size_factor <- SummarizedExomePeaks$sizeFactor

## Report the GLM statistics
SummarizedExomePeaks <- glmM(SummarizedExomePeaks) 

result =exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           gff_dir = GENE_ANNO_GTF,
           genome = "hg19",
           consistent_peak = TRUE,
           save_dir = direct_name,
           paired_end = FALSE)

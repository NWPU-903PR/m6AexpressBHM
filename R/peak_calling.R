peak_calling <- function(IP_BAM,
                         INPUT_BAM,
                        GENE_ANNO_GTF,
                        paired_end = FALSE,
                        Genome = "hg19"
                        output_dir){
##MeRIP-seq data alignment
MeRIP_Seq_Alignment <- scanMeripBAM(
                         bam_ip = IP_BAM,
                         bam_input = INPUT_BAM,
                         paired_end = FALSE)

##Perform peak calling on exons
###If either genome or bsgenome arguments are provided, the GC content bias correction will be performed while peak calling.
SummarizedExomePeaks <- exomePeakCalling(merip_bams = MeRIP_Seq_Alignment,
                                         gff_dir = GENE_ANNO_GTF,
                                         genome = Genome) 

##Estimate correction factors
###the normalization methods used will be controlled by the functions estimateSeqDepth() and normalizeGC(). 
####User can combine multiple normalization approaches, such as adding the conditional quantile normalization introduced in
SummarizedExomePeaks <- estimateSeqDepth(SummarizedExomePeaks)
SummarizedExomePeaks <- normalizeGC(SummarizedExomePeaks)
size_factor <- SummarizedExomePeaks$sizeFactor
## Report the GLM statistics
SummarizedExomePeaks <- glmM(SummarizedExomePeaks) 
## Export peak calling result
exportResults(SummarizedExomePeaks,
              save_dir = output_dir)
save(size_factor,file=paste0(output_dir,"/","size_factor.Rdata"))
}



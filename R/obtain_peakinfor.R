obtain_peakinfor <- function(peak_infor_dir){
  f1 <- paste0(peak_infor_dir,"/","Mod.csv")
  peak_sites_infor <- read.csv(f1)
  f2 <- paste0(peak_infor_dir,"/","ADDInfo","/","ADDInfo_ReadsCount.csv")
  peak_sites_reads <- read.csv(f2)
  seqnames <- as.character(peak_sites_infor$chr)
  start <- as.numeric(as.character(peak_sites_infor$chromStart))
  end <-  as.numeric(as.character(peak_sites_infor$chromEnd))
  strand <- as.character(peak_sites_infor$strand)
  gene_name <- as.character(peak_sites_infor$geneID)
  padj <- as.numeric(as.character(peak_sites_infor$padj))
  log2FoldChange <- as.numeric(as.character(peak_sites_infor$log2FoldChange))
  peak_infors <- data.frame(seqnames=seqnames,start=start,end=end,strand=strand,gene_name=gene_name,
                            peak_sites_reads[,-1],log2FoldChange=log2FoldChange,padj=padj)

  peaks_site_infors <- peak_infors
  return(peaks_site_infors)
}

peak_infors <- function(sites_infor,peak_sites_reads){
  seqnames <- as.character(peak_sites_infor$chr)
  start <- as.numeric(as.character(peak_sites_infor$chromStart))
  end <-  as.numeric(as.character(peak_sites_infor$chromEnd))
  strand <- as.character(peak_sites_infor$strand)
  gene_name <- as.character(peak_sites_infor$geneID)
  padj <- as.numeric(as.character(peak_sites_infor$padj))
  log2FoldChange <- as.numeric(as.character(peak_sites_infor$log2FoldChange))
  peak_infors <- data.frame(seqnames=seqnames,start=start,end=end,strand=strand,gene_name=gene_name,
                            peak_sites_reads[,-1],log2FoldChange=log2FoldChange,padj=padj)

  peaks_site_infor <- peak_infors
  return(peaks_site_infor)
}

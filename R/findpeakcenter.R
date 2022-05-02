###########obtain peak center
findpeakcenter <- function(annotation_file,maplongTX_peak){
  ##get the longest transcript
  txdbfile <- GenomicFeatures::makeTxDbFromGFF(annotation_file)
  exbytx_txdb <- exonsBy(txdbfile,by = "tx")
  isoform_ambiguity_method = "longest_tx"
  if(isoform_ambiguity_method == "longest_tx"){
    Longest_tx_table <- find_longest_transcript(exbytx_txdb,txdbfile)
    Kept_tx_indx <- Longest_tx_table$TXID[Longest_tx_table$longest]
    rm(Longest_tx_table)
  } else {
    Kept_tx_indx <- T
  }
  exbytx_txdb <- exbytx_txdb[Kept_tx_indx]
  exbytx_txdb <- exbytx_txdb[countOverlaps(exbytx_txdb,exbytx_txdb) == 1]
  ####
  targetpeaks <- maplongTX_peak$mapped_peankinfor
  maplongtx_peak <- maplongTX_peak$mapped_peakGRList
  #####
  targetpeak_GRlist <- .get_GRList(target_peak=targetpeaks,allpeak_GR=maplongtx_peak)
  targetpeak_center <- data.frame()
  txdbfiles <- exbytx_txdb
  for (i in 1:length(targetpeak_GRlist)) {
    onepeaksite <- unlist(targetpeak_GRlist[i])
    names(onepeaksite) <- NULL
    target_sites_map <-   mapToTranscripts(onepeaksite, txdbfiles,ignore.strand=F)
    if(length(target_sites_map)==1){
      peakcenter <- start(target_sites_map)+round((end(target_sites_map)-start(target_sites_map)+1)/2)
      peakcenterGR <- GRanges(seqnames = unique(names(txdbfiles[target_sites_map$transcriptsHits])),
                              IRanges(start = peakcenter,width = 1, names=unique(names(txdbfiles[target_sites_map$transcriptsHits]))),
                              strand = unique(as.character(strand(onepeaksite))))
      peakcenter_sites <- mapFromTranscripts(peakcenterGR,txdbfiles[target_sites_map$transcriptsHits])      
      onepeakcenter_infor <- data.frame(peaknum=as.character(names(targetpeak_GRlist[i])),seqnames=as.character(seqnames(peakcenter_sites)),
                                        start=start(peakcenter_sites),width=1,strand=as.character(strand(peakcenter_sites)),
                                        gene_name=unique(as.character(onepeaksite$gene_name)))
      targetpeak_center <- rbind(targetpeak_center,onepeakcenter_infor)
    }
    if(length(target_sites_map)>1){
      peakcenter <- start(target_sites_map)[1]+round((end(target_sites_map)[length(target_sites_map)]-start(target_sites_map)[1]+1)/2)
      peakcenterGR <- GRanges(seqnames = unique( names(txdbfiles[target_sites_map$transcriptsHits])),
                              IRanges(start = peakcenter,width = 1, names=unique(names(txdbfiles[target_sites_map$transcriptsHits]))),
                              strand = unique(as.character(strand(onepeaksite))))
      peakcenter_sites <- mapFromTranscripts(peakcenterGR,txdbfiles[unique(target_sites_map$transcriptsHits)])      
      onepeakcenter_infor <- data.frame(peaknum=as.character(names(targetpeak_GRlist[i])),seqnames=as.character(seqnames(peakcenter_sites)),
                                        start=start(peakcenter_sites),width=1,strand=as.character(strand(peakcenter_sites)),
                                        gene_name=unique(as.character(onepeaksite$gene_name)))
      targetpeak_center <- rbind(targetpeak_center,onepeakcenter_infor)
    }
    
  }
  return(targetpeak_center)
}
###
.get_GRList <- function(target_peak,allpeak_GR){
  target_peak_GR <- GRanges(seqnames = as.character(target_peak$seqnames),
                            IRanges(start = as.numeric(as.character(target_peak$start)),
                                    end = as.numeric(as.character(target_peak$end))),
                            strand = as.character(target_peak$strand),
                            gene_name=as.character(target_peak$gene_name))
  
  select_peaks <- data.frame()
  for (i in 1:length(target_peak_GR)) {
    one_overlap <- findOverlaps(allpeak_GR,target_peak_GR[i],type = "within")
    if(length(allpeak_GR[unique(one_overlap@from)])>1){
      one_selectpeak <- allpeak_GR[unique(one_overlap@from)[as.numeric(which.max(sum(width(allpeak_GR[unique(one_overlap@from)]))))]]
    }
    if(length(allpeak_GR[unique(one_overlap@from)])==1){
      one_selectpeak <- allpeak_GR[unique(one_overlap@from)]
    }
    one_select_peak <- unlist(one_selectpeak)
    one_select_peak$gene_name <- as.character(target_peak_GR$gene_name)[i]
    onepeakname <-names(one_select_peak)
    names(one_select_peak) <- NULL
    onepeak <- as.data.frame(one_select_peak)
    onepeak <- data.frame(peak_name=onepeakname,onepeak)
    select_peaks <- rbind(select_peaks,onepeak)
    
  }
  select_peakGR <- GRanges(seqnames = as.character(select_peaks$seqnames),
                           IRanges(start = as.numeric(as.character(select_peaks$start)),
                                   end = as.numeric(as.character(select_peaks$end))),
                           strand = as.character(select_peaks$strand))
  select_peakGR$exon_id <- as.numeric(as.character(select_peaks$exon_id))
  select_peakGR$exon_rank <- as.numeric(as.character(select_peaks$exon_rank))
  select_peakGR$gene_name <- as.character(select_peaks$gene_name)
  names(select_peakGR) <- as.character(select_peaks$peak_name)
  selectpeakGRlist <- split(select_peakGR,names(select_peakGR))
  return(selectpeakGRlist)
}

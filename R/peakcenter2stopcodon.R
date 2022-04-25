dist_stopcodon <- function(target_peakcenter,annotation_file){
  ##get the stop codon site
  txdbfile <- GenomicFeatures::makeTxDbFromGFF(annotation_file)
  genes_txdb <- genes(txdbfile)
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
  txid <- names(exbytx_txdb)
  cds <- cdsBy(txdbfile,by = "tx")
  cds <- cds[ names(cds) %in% txid ]
  Stop_codons <- resize( unlist( range(cds) ) , 3, fix = "end" )
  Stopcodon <- resize(Stop_codons,1,fix = "center")
  ##peak center
  target_sites_se <- GRanges(seqnames = as.character(target_peakcenter$seqnames),
                             IRanges(start = as.numeric(as.character(target_peakcenter$start)),
                                     width = 1),strand = as.character(target_peakcenter$strand))
  
  target_sites_map <-   mapToTranscripts(target_sites_se, exbytx_txdb,ignore.strand=F)
  stopcodon_map <- mapToTranscripts(Stopcodon,exbytx_txdb)
  ##
  overlap_hit_tx <- intersect(as.character(seqnames(target_sites_map)),
                              as.character(seqnames(stopcodon_map)))
  
  select_targetsites_map <- target_sites_map[which(!is.na(match(as.character(seqnames(target_sites_map)),overlap_hit_tx)))]
  select_peaknum <- as.numeric(as.character(select_targetsites_map$xHits))
  fun_hit <- function(overlap_hit){
    site_label <- which(!is.na(match(as.character(seqnames(select_targetsites_map)),overlap_hit)))
    stopcodon_label <- which(!is.na(match(as.character(seqnames(stopcodon_map)),overlap_hit)))
    site_pos <- start(select_targetsites_map[site_label])
    stopcodon_pos <- start(stopcodon_map[stopcodon_label])
    # tx_width <- as.numeric(sum(width(exbytx_txdb[which(!is.na(match(names(exbytx_txdb),overlap_hit)))])))
    abs_pos <- abs(site_pos-stopcodon_pos)
    # norm_pos <- abs_pos/tx_width
    site_num <- select_targetsites_map[site_label]$xHits
    names(abs_pos) <- site_num
    return(abs_pos)
  }
  
  site_norm_pos <-  mapply(fun_hit, overlap_hit_tx)
  names(site_norm_pos) <- NULL
  m6Asite_norm_pos <- unlist(site_norm_pos)
  m6Asite_num <- as.numeric(as.character(names(m6Asite_norm_pos)))
  select_m6Asite_overlap <- target_peakcenter[m6Asite_num,]
  add_norm_pos_ifor <- data.frame(select_m6Asite_overlap,as.numeric(as.character(m6Asite_norm_pos)))
  colnames(add_norm_pos_ifor)[ncol(add_norm_pos_ifor)] <- "dist_stopcodon"
  return(add_norm_pos_ifor)
}

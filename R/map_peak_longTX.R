##peak sites map to transcrpt

map_peak_longTX <- function(filepath,annotation_file,peak_sites_infor){
  ##obtain the longest transcript
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
  ##mapped to the longest transcript
  f1 <- paste0(filepath,"/","Mod.bed")
  a = read.table(f1,sep="\t",header=FALSE,stringsAsFactors =FALSE)
  read_peak <- import(filepath)
  mappeak_to_longtx <-  mapToTranscripts(read_peak, exbytx_txdb,ignore.strand=F)
  select_peak_label <- as.numeric(as.character(mappeak_to_longtx$xHits))
  a = a[select_peak_label,1:12]
  select_rawpeaks <- read_peak[select_peak_label]
  selectpeak_genomic <- mapFromTranscripts(mappeak_to_longtx,exbytx_txdb)
  select_peaks <- selectpeak_genomic[countOverlaps(selectpeak_genomic,selectpeak_genomic) == 1]
  last_select_label <- as.numeric(as.character(select_peaks$xHits))
  # mcols_info =a[,13:length(a[1,])]
  a = a[last_select_label,1:12]
  #############
  start_overlap_label <- which(!is.na(match(((a$V2)+1),peak_sites_infor$start)))
  a = a[start_overlap_label,]
  #################
  # get transcripts
  no_tx = length(a[,1])
  tx_id = 1:no_tx;
  tx_name = paste("peak",1:no_tx,sep="")
  tx_chrom = a[,1]
  tx_strand = a[,6]
  tx_start = a[,2]+1
  tx_end = a[,3]
  transcripts= data.frame(tx_id,tx_name,tx_chrom,tx_strand,tx_start,tx_end)
  
  # get genes
  tx_name = tx_name
  gene_id = as.character(a[,4])
  gene_id[is.na(gene_id)]="NA"
  gene=data.frame(tx_name,gene_id)
  
  # 
  splicing <- lapply(1:no_tx, .spliceSingleTrans, a=a, tx_start=tx_start)
  splicing <- .combineListOfSplicing(splicing)
  
  # make txdb
  peaks = suppressWarnings(
    makeTxDb(transcripts=transcripts, 
             splicings=splicing,
             genes=gene))
  
  # generate GRangesList
  tx <- exonsBy(peaks, "tx",use.names=TRUE)
  mcols(tx) <- a
  
  return(tx)
}

.combineListOfSplicing <- function(t){
  
  a <- paste("t[[",1:length(t),"]]", sep="")
  a <- paste(a,collapse =",")
  a <- paste("rbind(",a,")",sep="")
  c <- parse(text=a)
  b <- suppressWarnings(eval(c))
  
  return(b)
}

.spliceSingleTrans <- function(i,a,tx_start) {
  tx = a[i,]
  tx_id = i
  exon_rank=1:as.integer(tx[10])
  
  # get start
  temp = as.integer(strsplit(as.character(tx[12]), ",")[[1]]) + tx_start[i]
  exon_start=temp
  
  # get end
  temp = as.integer(strsplit(as.character(tx[11]), ",")[[1]])
  temp2 = temp + exon_start - 1
  exon_end=temp2
  
  # get CDS
  cds_start = exon_start
  cds_end = exon_end
  
  # get data frame
  splicing_tx = data.frame(tx_id,exon_rank,exon_start,exon_end,cds_start,cds_end)
  return(splicing_tx)
}

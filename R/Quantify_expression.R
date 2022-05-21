obtain_gene_expression <- function(Input_data, GTF_file, nthreads=20, isPairedEnd=F,Group_label,rep_num){
  gene_count <- featureCounts(Input_data,isGTFAnnotationFile=TRUE,GTF.featureType="exon",
                              GTF.attrType="gene_name", annot.ext = GTF_file, 
                              isPairedEnd=isPairedEnd, nthreads=nthreads)
 counts_data <- gene_count$counts
  countMatrix <- sapply(as.matrix(counts_data), as.numeric)
  countmatrixgene <- matrix(countMatrix, nrow=nrow(counts_data), ncol = ncol(counts_data))
  rownames(countmatrixgene) <- rownames(counts_data)
  sample_infor <- vector()
  for (i in 1:length(Group_label)) {
    one_infor <- paste0(Group_label[i],"_express",1:rep_num[i])
    sample_infor <- c(sample_infor,one_infor)
  }
  colnames(countmatrixgene) <- sample_infor
  
  condition <- rep(Group_label,rep_num)

  if(isPairedEnd==F){
    types=rep("single-read",length(condition))
  }
  if(isPairedEnd==T){
    types=rep("paired-read",length(condition))
  }
  coldata <- data.frame(condition=condition,type=types)
  coldata$condition <- factor(coldata$condition)
  coldata$type <- factor(coldata$type)
  rownames(coldata) <- sample_infor
  dds <- DESeqDataSetFromMatrix(countData = countmatrixgene,
                                colData = coldata,
                                design = ~ 1)
  size_factor <-  sizeFactors(estimateSizeFactors(dds)) 
  gene_expression_infor <- list(gene_express=countmatrixgene,
                                size_factor=size_factor)
  return(gene_expression_infor)   
}

  



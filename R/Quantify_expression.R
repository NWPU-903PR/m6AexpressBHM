obtain_gene_expression <- function(Input_data, GTF_file, nthreads=20, isPairedEnd=F){
  gene_count <- featureCounts(Input_data=Input_data,isGTFAnnotationFile=TRUE,GTF.featureType="exon",
                              GTF.attrType="gene_name", annot.ext = GTF_file, 
                              isPairedEnd=isPairedEnd, nthreads=nthreads)
counts_data <- gene_count$counts
countMatrix <- sapply(as.matrix(counts_data), as.numeric)
countmatrixgene <- matrix(countMatrix, nrow=nrow(counts_data), ncol = ncol(counts_data))
colnames(countmatrixgene) <- colnames(counts_data)
rownames(countmatrixgene) <- rownames(counts_data)
return(countmatrixgene)  
}

  



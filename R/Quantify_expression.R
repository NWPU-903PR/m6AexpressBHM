obtain_gene_expression <- function(Input_data, GTF_file, nthreads=20, isPairedEnd=F){
  gene_count <- featureCounts(Input_data=Input_data,isGTFAnnotationFile=TRUE,GTF.featureType="exon",
                              GTF.attrType="gene_name", annot.ext = GTF_file, 
                              isPairedEnd=isPairedEnd, nthreads=nthreads)
counts_data <- gene_count$counts
countMatrix <- sapply(as.matrix(counts_data), as.numeric)
countmatrixgene <- matrix(countMatrix, nrow=nrow(counts_data), ncol = ncol(counts_data))
colnames(countmatrixgene) <- colnames(counts_data)
rownames(countmatrixgene) <- rownames(counts_data)
conds <- factor(colnames(counts_data))
cds <- newCountDataSet(countmatrixgene, conds )
size_factor <-  sizeFactors(estimateSizeFactors(cds)) 

gene_expression_infor <- list(gene_express=countmatrixgene,
                              size_factor=size_factor)
return(gene_expression_infor)  
}

  



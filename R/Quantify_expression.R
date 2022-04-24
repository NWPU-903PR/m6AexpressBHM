gene_count <- featureCounts(Input_data,isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_name", annot.ext = gtf, isPairedEnd=T, nthreads=20)
counts_data <- gene_count$counts
colnames(counts_data) <- c("Im_rep1","Ma_rep1")
countMatrix <- sapply(as.matrix(counts_data), as.numeric)
countmatrixgene <- matrix(countMatrix, nrow=nrow(counts_data), ncol = ncol(counts_data))
colnames(countmatrixgene) <- colnames(counts_data)
rownames(countmatrixgene) <- rownames(counts_data)

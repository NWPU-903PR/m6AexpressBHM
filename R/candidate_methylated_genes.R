candidate_methylated_genes <- function(group_methylation,method){
Group_num <- length(group_methylation)
gene_names <- list()
for(i in 1:Group_num){
 gene_names[[i]] <- as.character(group_methylation[[i]]$gene_name)
}
consis_gene  <- Reduce(intersect, gene_names)
consis_gene_methy <- list()
for (i in 1:Group_num) {
  one_genemethy <- Group_methylation[[i]]
  onedata_consisgenemethy <- data.frame()
  for (j in 1:length(consis_gene)) {
    oneconsis_genemethy <- one_genemethy[one_genemethy$gene_name==consis_gene[j],-1]
    onedata_consisgenemethy <- rbind(onedata_consisgenemethy,oneconsis_genemethy)
  }
 colnames(onedata_consisgenemethy) <- paste0(colnames(onedata_consisgenemethy),"_Group",i)
 rownames(onedata_consisgenemethy) <- consis_gene
 consis_gene_methy[[i]] <- onedata_consisgenemethy
}
consis_genemethy <- do.call("cbind",consis_gene_methy)

}

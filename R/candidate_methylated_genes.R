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
if(method=="MAD"){
  within_group_methy <- data.frame()
  for (i in 1:Group_num) {
    Group_name <- paste0("Group",i)
    onegroup_methy <- consis_genemethy[,grep(Group_name,colnames(consis_genemethy ))]
    onecase_methy <- apply(onegroup_methy, 1, mean)
    within_group_methy <- rbind(within_group_methy,onecase_methy)
  }
  within_group_methy <- t(within_group_methy)
  rownames(within_group_methy) <- consis_gene
  colnames(within_group_methy) <- paste0("Group",1:Group_num,"_methy")
  methy_mean<- round(apply(within_group_methy, 1, mean),4) 
  methy_sd <- round(apply(within_group_methy, 1, sd),4)
  
  
}

}

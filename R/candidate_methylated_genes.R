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
  #methy_mean<- round(apply(within_group_methy, 1, mean),4) 
  #methy_sd <- round(apply(within_group_methy, 1, sd),4)
  methy_median <- apply(within_group_methy, 1, median)
  median_sd <- vector()
 for(i in 1:nrow(within_group_methy)){
   median_sd[i] <- median(abs(as.numeric(as.character(within_group_methy[i,])) - methy_median[i]))
 }
 MAD <- median_sd
 names(MAD) <- rownames(within_group_methy)
 select_genes <- names(MAD)[which(MAD>0.3)]
 candidate_gene_methy <- consis_genemethy[rownames(consis_genemethy)%in%select_genes,]
}
if(method=="DM"){
control_group <- consis_genemethy[,grep("Group1",colnames(consis_genemethy))]
design <- data.frame(Grp1=1,Grp2vs1=c(rep(0,ncol(control_group)*1),rep(1,(ncol(consis_genemethy)-ncol(control_group)))))
y <- consis_genemethy
# dupcor <- duplicateCorrelation(y,design,block=block)
fit1 <- lmFit(y,design)
fit1 <- eBayes(fit1)
diff_result <- topTable(fit1,coef = 2,number = nrow(y),genelist = as.character(rownames(consis_genemethy)))
DM_gene <- diff_result[diff_result$P.Value<0.05,]
DM_genename <- as.character(DM_gene$ID)
candidate_gene_methy <-  consis_genemethy[rownames(consis_genemethy)%in%DM_genename,]
}
return(candidate_gene_methy)
}

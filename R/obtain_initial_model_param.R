initial_parameters <- function(input_data_path,input_data_name,output_path,output_name,num_group,size_factor){
  fa <- paste0(input_data_path,input_data_name)
  gene_count_methy <- read.table(fa,header = T)
  
  gene_reads <- gene_count_methy[,2:((ncol(gene_count_methy)+1)/2)]
  select_gene <- which(apply(gene_reads, 1, mean)>10)
  select_counts_methy <- gene_count_methy[select_gene,]
  gene_name <- as.character(select_counts_methy$gene_name)
  # select_counts_methy <- gene_count_methy
  cl<- makeCluster(detectCores()-28)      
  registerDoParallel(cl) 
  suppressMessages(suppressWarnings(glmmtmb_results<-foreach(i=1:nrow(select_counts_methy),.combine = rbind,.packages = "glmmTMB") %dopar% {
    one_genecount <- as.numeric(select_counts_methy[i,2:((ncol(select_counts_methy)+1)/2)])
    one_genemethy <- as.numeric(select_counts_methy[i,(((ncol(select_counts_methy)+1)/2)+1):ncol(select_counts_methy)])
    y=one_genecount
    sample_ID <- seq(1,length(one_genecount),1)
    subject_ID <- rep(seq(1,num_group,1),rep(length(one_genecount)/num_group,num_group))
    factors <- data.frame(subject_ID,sample_ID,one_genemethy)
    subject =(factors[, "subject_ID"])
    methy <- factors$one_genemethy
    suppressMessages(res <- try(glmmTMB(y~ methy+(1|subject),family = nbinom2,data = factors,offset = size_factor)))
    if(inherits(res, "try-error")){
      next
    }else{
      suppressMessages(suppressWarnings(f0 <- glmmTMB(y~ methy+(1|subject),family = nbinom2,data = factors,offset = size_factor)))
      model_summary <- summary(f0)
      fixed_coefinfor <- model_summary[["coefficients"]][["cond"]]
      beta <- as.numeric(as.character(fixed_coefinfor[,1]))
      pvalue <- as.numeric(as.character(fixed_coefinfor[,4]))
      alpha_value <- 1/model_summary[["sigma"]]
      one_coef_infor <- c(beta,pvalue,alpha_value)
      b_variance <- as.numeric(model_summary[["varcor"]][["cond"]][["subject"]])
      parfull_data <- f0[["fit"]][["parfull"]]
      one_b <- c(parfull_data[which(names(parfull_data)=="b")],b_variance)
      names(one_b) <- NULL
      one_data <- c((gene_name[i]),one_coef_infor,one_b)
    }
  }
  ))
  stopCluster(cl )
  colnames(glmmtmb_results) <- c("gene_name","beta0","beta1","beta0_pvalue","beta1_pvalue","alpha",
                                 paste0("b0_",1:num_group),"variance_b0")
  rownames(glmmtmb_results) <- NULL
  glmmtmb_results <- as.data.frame(glmmtmb_results)
  writexl::write_xlsx(glmmtmb_results,paste0(output_path,"glmmtmb_result_peaklevel.xlsx"))
  # return(glmmtmb_results)
  select_genelabel <- glmmtmb_results$gene_name
  select_genes <-as.character(select_genelabel)
  coef_infor <- glmmtmb_results[,2:6]
  b_infor <- glmmtmb_results[,7:ncol(glmmtmb_results)]
  select_label <- union(which(!is.na(coef_infor$beta0_pvalue)),which(!is.na(coef_infor$beta1_pvalue)))
  select_b_infor <- b_infor[select_label,]
  select_coefinfor <- coef_infor[select_label,]
  select_gene_name <- (select_genes[select_label])
  output_data <- list(select_b_infor,select_coefinfor,select_gene_name)
  save(output_data,file = paste0(output_path, output_name))
}

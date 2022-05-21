##Match gene expression and methylation level for candidate genes
match_expr_methy <- function(gene_express_data,gene_methy){
gene_express <- gene_express_data$gene_express
express_sizefactor <- gene_express_data$size_factor
full_expression <- apply(gene_express,1,sum)
select_label <- which(full_expression>10)
full_expr_gene <- gene_express[select_label,]
target_genes <- intersect(rownames(gene_methy),rownames(full_expr_gene))
gene_exp_methy <- data.frame()
for (i in 1:length(target_genes)) {
  one_exp_methy <- cbind(full_expr_gene[rownames(full_expr_gene)==target_genes[i],],
                         gene_methy[rownames(gene_methy)==target_genes[i],])
  gene_exp_methy <- rbind(gene_exp_methy,one_exp_methy)
}
 gene_expr_methy_infor <- list(gene_expr_methy=gene_exp_methy,
                               library_sizefactor=express_sizefactor)
 return(gene_expr_methy_infor)
}

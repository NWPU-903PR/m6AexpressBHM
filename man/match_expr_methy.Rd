\name{match_expr_methy}
\alias{match_expr_methy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Match expression and methylation
}
\description{
Match the paired expression and methylation level for each gene.
}
\usage{
match_expr_methy(gene_express_data,gene_methy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_express_data}{
a \code{list} containing a \code{matrix} with gene expression for each sample in each group and a \code{numeric} vector specifying the sequencing depth size factors for each Input samples.
  \item{gene_methy}{
a \code{data.frame} containing the candidate genes' methylation for each sample in each group.
}
}

\value{
It will return a list containing two components:
\item{gene_expr_methy}{a \code{data.frame} including the matched gene expression and methylation level for each sample in all concerned groups.
\item{size_factor}{a \code{numeric} vector specifying the sequencing depth size factors for each Input samples of all concerned groups.}
}
\dontrun{
\examples{
### Match expression and methylation genes
gene_expr_methy <- match_expr_methy(gene_express_data=gene_expression,                                                     gene_methy=select_methylated_genes)
}
}

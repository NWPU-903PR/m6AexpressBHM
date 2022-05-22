\name{candidate_methylated_genes}
\alias{candidate_methylated_genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select candidate methylation genes
}
\description{
Select candidate methylation genes from multiple concerned groups.
}
\usage{
candidate_methylated_genes(group_methylation,
                           method,
                           threshod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group_methylation}{
A \code{list} includes the methylation level of each genes obtained by \code{\link{gene_methy_level_distdecay}} function from multiple concerned groups.
}
  \item{method}{
A \code{character} specifies the method of selecting candidate genes, which incudes median absolute deviation (\code{"MAD"}) method, differential methylation method (\code{"DM"}) and obtaines longitudinal methylation genes (\code{"longitudinal methylation"}).
}
  \item{threshod}{
A \code{numeric} specifies the threshod for different methods to select candidate genes. In \code{method="MAD"}, the \code{threshod} is set \code{"threshod=0.3"} by default. In \code{method="DM"} context, the \code{threshod} is the \code{p-value} threshod setting as \code{"threshod=0.05"} by default. In \code{"longitudinal methylation"} method, the \code{threshod} is the \code{p-value} threshod setting as \code{"threshod=0.05"} by default to select time correlated methylation genes.
}
}

\value{
It will report a \code{data.frame} including the methylation level of each candidate genes in multiple groups.
}

\dontrun{
\examples{
###select candidate methylated genes with high variable methylation genes (e.g. MAD>0.3)
Group_methylation <- list(Group1_methy=Group1_methylevel_distdecay,
                          Group2_methy=Group2_methylevel_distdecay,
                          Group3_methy=Group3_methylevel_distdecay)
select_methylated_genes <- candidate_methylated_genes(group_methylation=Group_methylation,
                                                      method="MAD",threshod=0.3)

}

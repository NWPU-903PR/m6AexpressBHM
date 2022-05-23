\name{obtain_initial_parameters}
\alias{obtain_initial_parameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain initial parameters
}
\description{
Obtain initial model parameters.
}
\usage{
obtain_initial_parameters(gene_expre_methy,
                          num_group,
                          size_factor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_expre_methy}{
a \code{data.frame} including the matched gene expression and methylation level for each sample in all concerned groups.
}
  \item{num_group}{
  a \code{integer} giving the number of concerned groups.
  }
  \item{size_factor}{
  a \code{numeric} vector specifying the sequencing depth size factors for each Input samples.
  }
}

\value{
It will return a \code{list} including three components:
  \item{b_infor}{a \code{data.frame} including the random effects value and their corresponding variance, which reflects the variation of gene expression across different concerned groups}
  \item{coeff_infor}{a \code{data.frame} inclduing the columns of the initially predicting \code{Beta0}, \code{Beta1}, \code{alpha}. \code{Beta0} is a gene-specific intercept that models the baseline log gene expression. \code{Beta1} captures the influece of m6A methylation on gene expression in all concerned groups. \code{alpha} is a gene-specific dispersion parameter in negative binomial regression model.}
  \item{gene_name}{a \code{character} vector inclduing the gene symbol or gene name}
}

\examples{
\dontrun{
##the match gene reads count and methylation level
gene_reads_methy <- gene_expr_methy$gene_expr_methy
##library size factor for each gene expression data
size_factor <- gene_expr_methy$library_sizefactor
##Obtain the initail model parameters
initial_model_param <- obtain_initial_parameters(gene_expre_methy=gene_reads_methy,
                                                 num_group=3,size_factor=size_factor)
}
}

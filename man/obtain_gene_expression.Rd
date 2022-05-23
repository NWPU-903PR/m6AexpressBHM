\name{obtain_gene_expression}
\alias{obtain_gene_expression}
\title{
Obtain gene expression
}
\description{
Obtain the reads count for each gene to represent the expression of each gene and report the sequencing depth size factors for each Input samples
}
\usage{
obtain_gene_expression(Input_data, 
                       GTF_file, 
                       nthreads=1, 
                       isPairedEnd=FALSE,
                       Group_label,
                       rep_num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Input_data}{
a \code{character} vector for the BAM file directories of the Input samples, which are used to quantify gene expression.
}
  \item{GTF_file}{
a \code{character} specifies the directory toward a gene annotation GFF/GTF file.
}
  \item{nthreads}{
a \code{integer} giving the number of threads used for running this function. 1 by default.
}
  \item{isPairedEnd}{
  a \code{logical} scalar or a \code{logical} vector, indicating whether libraries contain paired-end reads or not. FALSE by default.
  }
  \item{Group_label}{
a \code{character} vector giving the name of each concerned group.
}
  \item{rep_num}{
a \code{integer} vector giving the number samples of each concerned group.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
It will return a list containing two components:
\item{gene_express}{a \code{matrix} including gene reads count in each sample of each concerned group}
\item{size_factor}{a \code{numeric} vector specifying the sequencing depth size factors for each Input samples.}
}

\examples{
\dontrun{
### obtain gene reads count from Input dataset in BAM format
Input_data <- c(group1_INPUT_BAM,group2_INPUT_BAM,group3_INPUT_BAM)
gene_expression <- obtain_gene_expression(Input_data=Input_data, 
                                          GTF_file=GENE_ANNO_GTF, 
                                          nthreads=10, 
                                          isPairedEnd=FALSE,
                                          Group_label=c("Group1","Group2","Group3"),
                                          rep_num=c(3,3,3))
}
}

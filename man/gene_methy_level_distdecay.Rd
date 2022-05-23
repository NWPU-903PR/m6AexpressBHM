\name{gene_methy_level_distdecay}
\alias{gene_methy_level_distdecay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain gene's methylation level
}
\description{
Calculate methylation level of each gene, which is weighted by the distance of peak center to stop-codon
}
\usage{
gene_methy_level_distdecay(mapLTX_peakinfor,
                           size_factor,
                           peak_dist_stopcodon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mapLTX_peakinfor}{
A \code{list} containing peak sites information mapped to the longest transcript of gene.
}
  \item{size_factor}{
A \code{numberic vector} containing the sequencing depth size factors for each MeRIP-seq samples, which are save in the \code{\link{peak_calling}} output directory in a \strong{"size_factor.Rdata"} file
}
  \item{peak_dist_stopcodon}{
A \code{data.frame} specifing the distance of peak center to stop cond and the site of each peak center.
}
}

\value{
It will return a \code{data.frame} file containing the methylation level weighted by the distance of peak center to sop codon for each gene in each concerned sample.
}

\examples{
\dontrun{
load("./Group1_peakcalling/size_factor.Rdata")
Group1_sizefactor <- size_factor
Group1_methylevel_distdecay <- gene_methy_level_distdecay(mapLTX_peakinfor=Group1_mappeak_LTX,
                                                          peak_dist_stopcodon=Group1_peakcenter2stopcondon,
                                                          size_factor=Group1_sizefactor)

}
}

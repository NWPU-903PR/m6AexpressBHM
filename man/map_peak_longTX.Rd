\name{map_peak_longTX}
\alias{map_peak_longTX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Map peak sites to the longest transcript
}
\description{
Map m6A peak sites to the longest transcirpt for each gene
}
\usage{
map_peak_longTX(filepath,
                annotation_file,
                peak_sites_infor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{
A \code{character} specifies the directory of the peak calling result obtained by \code{\link{peak_calling}} function.
}
  \item{annotation_file}{
A \code{character} specifies the directory toward a gene annotation GFF/GTF file.
}
  \item{peak_sites_infor}{
A \code{data.frame} file gives the m6A peak sites information obtained by \code{\link{obtain_peakinfor}} function.
}
}
\value{
A list containing two components:
  \item{mapped_peakGRList}{peak sites mapped to the longest transicript in \code{GRangesList} format}
  \item{mapped_peankinfor }{peak sites information mapped to the longest transcirpt of each gene in \code{data.frame} format }
}
\dontrun{
\examples{
Group1_mappeak_LTX <- map_peak_longTX(filepath="./Group1_peakcalling",
                                     annotation_file=GENE_ANNO_GTF,
                                     peak_sites_infor=Group1_peakinfor)
}

\name{findpeakcenter}
\alias{findpeakcenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain peak center
}
\description{
Obtain each peak center site based on transcriptome coordinate.
}
\usage{
findpeakcenter(annotation_file,
              maplongTX_peak)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotation_file}{
A \code{character} specifies the directory toward a gene annotation GFF/GTF file.
}
  \item{maplongTX_peak}{
A \code{list} includes the information of peak sites mapped to the longest transcript of each gene.
}
}
\value{
It will return a \code{data.frame} file including the sites of peak center for each peak.
}

\dontrun{\examples{
Group1_peakcenter <- findpeakcenter(annotation_file=GENE_ANNO_GTF,
                                    maplongTX_peak=Group1_mappeak_LTX)
}

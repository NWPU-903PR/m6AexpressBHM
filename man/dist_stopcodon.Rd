\name{dist_stopcodon}
\alias{dist_stopcodon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance to stop codon
}
\description{
Calculate the distance between peak center to stop condon.
}
\usage{
dist_stopcodon(target_peakcenter,
              annotation_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target_peakcenter}{
A \code{data.frame} file including the site of peak center for each peak.
}
  \item{annotation_file}{
A \code{character} specifies the directory toward a gene annotation GFF/GTF file.
}
}
\details{
The calculation of distance between the peak center and stop codon is based on the transcriptome coordinates eliminating the disturbance of intron.
}
\value{
It will return a \code{data.frame} file specifing the distance of peak center to stop cond and the site of each peak center.
}

\dontrun{
\examples{
Group1_peakcenter2stopcondon <- dist_stopcodon(target_peakcenter=Group1_peakcenter,
                                                annotation_file=GENE_ANNO_GTF)
}
}

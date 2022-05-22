\name{peak_calling}
\alias{peak_calling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Peak calling based on \code{\link{exomePeak2}} for MeRIP-seq dataset
}
\usage{
peak_calling(IP_BAM = NULL,
             INPUT_BAM = NULL,
             GENE_ANNO_GTF = NULL,
             paired_end = FALSE,
             Genome = NA,
             output_dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IP_BAM}{a \code{character} vector for the BAM file directories of the IP samples.
}
  \item{INPUT_BAM}{a \code{character} vector for the BAM file directories of the Input samples.
}
  \item{GENE_ANNO_GTF}{a \code{character} specifies the directory toward a gene annotation GFF/GTF file.
}
  \item{Genome}{a \code{\link{BSgenome}} object for the genome sequence information. For example: \code{Genome = "hg19"}.}
  \item{output_dir}{a \code{character} for the name of the directory being output; Default \code{= "exomePeak2_output"}.}
}
\description{
\code{peak_calling} conducts peak calling and peak statistics calculation from \strong{BAM} files of a MeRIP-seq experiment.
\enumerate{
\item Check and index the BAM files with \code{\link{scanMeripBAM}}.
\item Call modification peaks on exons with \code{\link{exomePeakCalling}}.
\item Estimate sequencing depth size factors for each MeRIP-seq samples  \code{\link{exomePeakCalling}}.
\item Calculate offset factors of GC content biases with \code{\link{estimateSeqDepth}}.
\item Calculate (differential) modification statistics with the generalized linear model (GLM) using \code{\link{glmM}}.
\item Export the peaks/sites statistics with user defined format by \code{\link{exportResults}}.
}
}
\details{
\code{\link{peak_calling}} call RNA modification peaks and calculate peak statistics from \strong{BAM} files of a MeRIP-seq experiment based on \code{\link{exomePeak2}} R package.

The transcript annotation from the \strong{GFF} file) should be provided to perform analysis on exons.

The \code{\link{BSgenome}} object is also required to perform the GC content bias adjustment.
If the \code{Genome} argument is not provided (\code{= NULL}), the downstream analysis will proceed without GC content bias corrections.

Under the default setting, \code{\link{peak_calling}} will save the results of modification analysis under a folder named \code{'exomePeak2_output'}.
The results generated include a \strong{BED} file and a \strong{CSV} table that stores the locations and statistics of the modified peaks/sites.
\code{peak_calling} also report the sequencing depth size factors for each MeRIP-seq samples, which are save in the same output directory in a \strong{"size_factor.Rdata"} file

}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 See Also as \code{\link{?exomePeak2}}
}
\dontrun{
\examples{
library(exomePeak2)
f1 <- "./group1_IP1.bam"
f2 <- "./group1_IP2.bam"
f3 <- "./group1_IP3.bam"
f4 <- "./group1_Input1.bam"
f5 <- "./group1_Input2.bam"
f6 <- "./group1_Input3.bam"
group1_IP_BAM <- c(f1,f2,f3)
group1_INPUT_BAM <- c(f4,f5,f6)
GENE_ANNO_GTF = "./hg19_GTF/genes.gtf"
###peak calling for group1
group1_peak_calling <- peak_calling(IP_BAM=group1_IP_BAM,
                                    INPUT_BAM=group1_INPUT_BAM,
                                    GENE_ANNO_GTF=GENE_ANNO_GTF,
                                    paired_end = FALSE,
                                    Genome = "hg19"
                                    output_dir="./Group1_peakcalling")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

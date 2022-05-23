\name{obtain_peakinfor}
\alias{obtain_peakinfor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain peak sites information.
}
\description{
Obtain peak sites information from the output of \code{\link{peak_calling}}.
}
\usage{
obtain_peakinfor(peak_infor_dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peak_infor_dir}{
a \code{character} specifies the directory of peak calling results by \code{\link{peak_calling}} function.
}
}
\details{
This function will output the m6A peak sites informaiton in a \code{data.frame} file, which includes the \code{seqnames}, \code{start}, \code{end}, \code{strand}, the \code{reads count} in all concerned samples of each peak site, the \code{log2Foldchagne} and \code{padj} for each peak.
}
\value{

  \item{peaks_site_infors}{A \code{data.frame} file, which includes the \code{seqnames}, \code{start}, \code{end}, \code{strand}, the \code{reads count} in all concerned samples of each peak site, the \code{log2Foldchagne} and \code{padj} for each peak.}


}


\examples{
\dontrun{
Group1_peakinfor <- obtain_peakinfor(peak_infor_dir="./Group1_peakcalling")
}
}

\name{peaks.gr}
\alias{peaks.gr}
\docType{data}
\title{
example cleavage sites
}
\description{
An example data set containing cleavage sites (peaks) from getPeaks
}
\usage{data("peaks.gr")}
\format{
    \describe{GRanges with count (peak height), bg (local background), 
        SNratio (signal noise ratio), p-value, and option adjusted p-value }
}
\value{
   \item{peaks.gr}{GRanges with count (peak height), bg (local background),
       SNratio (signal noise ratio), p-value, and option adjusted p-value }
}
\source{
http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1695644
}
\examples{
    data(peaks.gr)
    names(peaks.gr)
    peaks.gr
}
\keyword{datasets}

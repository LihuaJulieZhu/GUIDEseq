\name{getPeaks}
\alias{getPeaks}
\title{
Obtain peaks from GUIDE-seq
}
\description{
Obtain strand-specific peaks from GUIDE-seq
}
\usage{
getPeaks(gr, window.size = 20L, step = 20L, bg.window.size = 5000L,
    min.reads = 10L, min.SNratio = 2, maxP = 0.05, 
    stats = c("poisson", "nbinom"), p.adjust.methods = 
    c("none", "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr"))
}
\arguments{
  \item{gr}{
GRanges with cleavage sites, output from getUniqueCleavageEvents
}
  \item{window.size}{
window size to calculate coverage
}
  \item{step}{
step size to calculate coverage
}
  \item{bg.window.size}{
window size to calculate local background
}
  \item{min.reads}{
minimum number of reads to be considered as a peak
}
  \item{min.SNratio}{
minimum signal noise ratio, which is the coverage normalized by local 
    background  
}
  \item{maxP}{
Maximum p-value to be considered as significant
}
  \item{stats}{
Statistical test, default poisson
}
  \item{p.adjust.methods}{
Adjustment method for multiple comparisons, default none
}
}
\value{
   \item{peaks }{GRanges with count (peak height), bg (local background), 
       SNratio (signal noise ratio), p-value, and option adjusted p-value }
   \item{summarized.count}{A data frame contains the same information
       as peaks except that it has all the sites without filtering.
    }
}
\author{
Lihua Julie Zhu
}

\examples{
    if (interactive())
    {
        data(uniqueCleavageEvents)
        peaks <- getPeaks(uniqueCleavageEvents$cleavage.gr,
            min.reads = 80, multicore = TRUE, n.cores.max = 2)
        peaks$peaks
    }
}
\keyword{misc}

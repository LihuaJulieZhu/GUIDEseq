\name{mergePlusMinusPeaks}
\alias{mergePlusMinusPeaks}
\title{
Merge peaks from plus strand and minus strand
}
\description{
Merge peaks from plus strand and minus strand with required orientation and 
    within certain distance apart}
\usage{
mergePlusMinusPeaks(peaks.gr, peak.height.mcol = "count", 
    bg.height.mcol = "bg", distance.threshold = 40L, 
    max.overlap.plusSig.minusSig = 10L,
    plus.strand.start.gt.minus.strand.end = TRUE, output.bedfile)
}
\arguments{
  \item{peaks.gr}{
Specify the peaks as GRanges object, which should contain peaks
    from both plus and minus strand. In addition, it should contain peak height
    metadata column to store peak height and optionally background height.
}
  \item{peak.height.mcol}{
Specify the metadata column containing the peak height, default to count
}
  \item{bg.height.mcol}{
Specify the metadata column containing the background height, default to bg
}
  \item{distance.threshold}{
Specify the maximum gap allowed between the plus stranded and 
    the negative stranded peak, default 40. Suggest set it to twice of
    window.size used for peak calling.
}
\item{max.overlap.plusSig.minusSig}{
Specify the cushion distance (maximum overlap between plus strand peak and minus strand peak).
Default to 10L to allow sequence error and inprecise integration. Only applicable
if plus.strand.start.gt.minus.strand.end is set to TRUE.
}
  \item{plus.strand.start.gt.minus.strand.end}{
Specify whether plus strand peak start greater than
    the paired negative strand peak end. Default to TRUE
}
  \item{output.bedfile}{
Specify the bed output file name, which is used for off target analysis 
    subsequently.
}
}
\value{
output a list and a bed file containing the merged peaks
    a data frame of the bed format
   \item{mergedPeaks.gr}{merged peaks as GRanges}
   \item{mergedPeaks.bed}{merged peaks in bed format}
}
\references{ 
Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to 
    annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 
    11:237doi:10.1186/1471-2105-11-237.
Zhu L.J. (2013) Integrative analysis of ChIP-chip and ChIP-seq dataset. 
    Methods Mol Biol. 2013;1067:105-24. doi: 10.1007/978-1-62703-607-8\_8.
}
\author{
Lihua Julie Zhu
}
\examples{

if (interactive())
{
    data(peaks.gr)
    mergedPeaks <- mergePlusMinusPeaks(peaks.gr = peaks.gr,
        output.bedfile = "mergedPeaks.bed")
    mergedPeaks$mergedPeaks.gr
    head(mergedPeaks$mergedPeaks.bed)
}
}
\keyword{ misc }

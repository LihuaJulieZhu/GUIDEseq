\name{getUsedBarcodes}
\alias{getUsedBarcodes}
\title{
Create barcodes from the p5 and p7 index used for each sequencing lane
}
\description{
Create barcodes from the p5 and p7 index for assigning reads to each barcode
}
\usage{
getUsedBarcodes(p5.index, p7.index, header = FALSE, reverse.p7 = TRUE, 
    reverse.p5 = FALSE, outputFile)
}
\arguments{
  \item{p5.index}{
  A text file with one p5 index sequence per line
}
  \item{p7.index}{
  A text file with one p7 index sequence per line
}
  \item{header}{
  Indicate whether there is a header in the p5.index and p7.index files.
  Default to FALSE
}
  \item{reverse.p7}{
  Indicate whether to reverse p7 index, default to TRUE for standard
GUIDE-seq experiments
}
  \item{reverse.p5}{
  Indicate whether to reverse p5 index, default to FALSE for standard 
GUIDE-seq experiments
}
  \item{outputFile}{
  Give a name to the output file where the generated barcodes are written
}
}
\value{
DNAStringSet
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lihua Julie Zhu
}
\note{
Create barcode file to be used to bin the reads sequenced in a pooled lane
}

\examples{
    p7 <- system.file("extdata", "p7.index",
           package = "GUIDEseq") 
    p5 <- system.file("extdata", "p5.index",
           package = "GUIDEseq")
    outputFile <- "usedBarcode" 
    getUsedBarcodes(p5.index = p5, p7.index = p7, reverse.p7 = TRUE,
        reverse.p5 = FALSE, outputFile = outputFile)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ utilities }

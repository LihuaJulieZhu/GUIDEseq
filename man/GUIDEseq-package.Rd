\name{GUIDEseq-package}
\alias{GUIDEseq-package}
\alias{GUIDEseq}
\docType{package}
\title{
Analysis of GUIDE-seq}
\description{
The package includes functions to retain one read per unique molecular 
    identifier (UMI), filter reads lacking integration oligo sequence, identify 
peak locations (cleavage sites) and heights, merge peaks, perform off-target 
search using the input gRNA. This package leverages CRISPRseek and ChIPpeakAnno
packages.
}
\details{
\tabular{ll}{
Package: \tab GUIDEseq\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-09-04\cr
License: \tab GPL (>= 2)\cr
}
Function GUIDEseqAnalysis integrates all steps of GUIDE-seq analysis into one 
function call
}
\author{
Lihua Julie Zhu
Maintainer:julie.zhu@umassmed.edu
}
\references{
Shengdar Q Tsai and J Keith Joung et al. GUIDE-seq enables genome-wide 
profiling of off-target cleavage by CRISPR-Cas nucleases. Nature Biotechnology
33, 187 to 197 (2015)
}
\keyword{ package }
\seealso{
GUIDEseqAnalysis
}
\examples{
    if(interactive())
    {
        library("BSgenome.Hsapiens.UCSC.hg19")
        umiFile <- system.file("extdata", "UMI-HEK293_site4_chr13.txt",
           package = "GUIDEseq")
        alignFile <- system.file("extdata","bowtie2.HEK293_site4_chr13.sort.bam" ,
            package = "GUIDEseq")
        gRNA.file <- system.file("extdata","gRNA.fa", package = "GUIDEseq")
        guideSeqRes <- GUIDEseqAnalysis(
            alignment.inputfile = alignFile,
            umi.inputfile = umiFile, gRNA.file = gRNA.file,
            orderOfftargetsBy = "peak_score",
            descending = TRUE,
            keepTopOfftargetsBy = "predicted_cleavage_score",
            scoring.method = "CFDscore",
            BSgenomeName = Hsapiens, min.reads = 80, n.cores.max = 1)
        guideSeqRes$offTargets
   }
}

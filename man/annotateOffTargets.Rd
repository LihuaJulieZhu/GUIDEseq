\name{annotateOffTargets}
\alias{annotateOffTargets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotate offtargets with gene name
}
\description{
Annotate offtargets with gene name and whether it is inside an exon
}
\usage{
annotateOffTargets(thePeaks, txdb, orgAnn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{thePeaks}{
Output from offTargetAnalysisOfPeakRegions 
}
\item{txdb}{
TxDb object, for creating and using TxDb object, please refer 
to GenomicFeatures package. For a list of existing TxDb object, please 
search for annotation package starting with Txdb at 
http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData, 
such as TxDb.Rnorvegicus.UCSC.rn5.refGene for rat, 
TxDb.Mmusculus.UCSC.mm10.knownGene for mouse, 
TxDb.Hsapiens.UCSC.hg19.knownGene for human, 
TxDb.Dmelanogaster.UCSC.dm3.ensGene for Drosophila and 
TxDb.Celegans.UCSC.ce6.ensGene for C.elegans
}
\item{orgAnn}{
organism annotation mapping such as org.Hs.egSYMBOL in org.Hs.eg.db package for human
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame and a tab-delimited file offTargetsInPeakRegions.xls, 
    containing all input offtargets
    with potential gRNA binding sites, mismatch number and positions, alignment
    to the input gRNA and predicted cleavage score, and whether the offtargets are
    inside an exon and associated gene name.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lihua Julie Zhu
}
\seealso{
GUIDEseqAnalysis
}
\examples{
if (!interactive()) {
    library("BSgenome.Hsapiens.UCSC.hg19")
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    peaks <- system.file("extdata", "T2plus100OffTargets.bed",
        package = "CRISPRseek")
    gRNAs <- system.file("extdata", "T2.fa",
        package = "CRISPRseek")
    outputDir = getwd()
    offTargets <- offTargetAnalysisOfPeakRegions(gRNA = gRNAs, peaks = peaks,
        format=c("fasta", "bed"),
        peaks.withHeader = TRUE, BSgenomeName = Hsapiens,
        upstream = 20L, downstream = 20L, PAM.size = 3L, gRNA.size = 20L,
        orderOfftargetsBy = "predicted_cleavage_score", 
        PAM = "NGG", PAM.pattern = "(NGG|NAG|NGA)$", max.mismatch = 2L,
        outputDir = outputDir,
        allowed.mismatch.PAM = 3, overwrite = TRUE)
    annotatedOfftargets <- annotateOffTargets(offTargets, 
       txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
       orgAnn = org.Hs.egSYMBOL)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}% use one of  RShowDoc("KEYWORDS")

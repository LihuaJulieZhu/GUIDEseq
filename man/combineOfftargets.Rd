\name{combineOfftargets}
\alias{combineOfftargets}
\title{
Combine Offtargets
}
\description{
Merge offtargets from different samples
}
\usage{
combineOfftargets(offtarget.folder, sample.name, 
    remove.common.offtargets = FALSE,control.sample.name,
    offtarget.filename = "offTargetsInPeakRegions.xls", 
    common.col = c("targetSeqName", "chromosome", 
        "offTargetStrand", "offTarget_Start", 
        "offTarget_End", "gRNAPlusPAM", "offTarget_sequence",
         "n.mismatch", "guideAlignment2OffTarget", 
         "predicted_cleavage_score"), 
    exclude.col = "name", 
    outputFileName)
}
\arguments{
  \item{offtarget.folder}{
   offtarget summary output folders created in GUIDEseqAnalysis function
}
  \item{sample.name}{
   Sample names to be used as part of the column names
   in the final output file
}
\item{remove.common.offtargets}{
If set to TRUE, off-targets common to all samples will be removed.
}
\item{control.sample.name}{
 The name of the control sample for filtering off-targets present in the control sample
}
  \item{offtarget.filename}{
   Default to offTargetsInPeakRegions.xls, generated in GUIDEseqAnalysis
      function
}
  \item{common.col}{
   common column names used for merge files. Default to 
       c("targetSeqName", "chromosome",
       "offTargetStrand", "offTarget_Start",
       "offTarget_End", "gRNAPlusPAM", "offTarget_sequence",
       "n.mismatch", "guideAlignment2OffTarget",
       "predicted_cleavage_score")
}
  \item{exclude.col}{
  columns to be excluded before merging. Default to name (second column of 
  offTargetsInPeakRegions.xls). Please check offTargetsInPeakRegions.xls
  to choose the desired columns to exclude
}
  \item{outputFileName}{
  The merged offtarget file
}
}
\details{
  Please note that by default, merged file will only contain peaks with offtargets
  found in the genome in GUIDEseqAnalysis function. 
}
\value{
    a tab-delimited file similar to offTargetsInPeakRegions.tsv, containing all peaks
    from all samples merged by potential gRNA binding sites, mismatch number and positions,
    alignment to the input gRNA and predicted cleavage score. Sample specific columns have
    sample.name concatenated to the original column name, e.g., peak_score becomes 
    sample1.peak_score.
}
\author{
Lihua Julie Zhu
}

\examples{
    offtarget.folder <- system.file("extdata", 
        c("sample2-18", "sample3-19", "sample4-20"),
        package = "GUIDEseq")
    mergedOfftargets <- 
       combineOfftargets(offtarget.folder = offtarget.folder, 
       sample.name = c("Wild-type SpCas9", "SpCas9-MT3-ZFP", "Split-SpCas9 dual NLS"), 
       outputFileName = "TS2offtargets3Constructs.xls")
}
\keyword{ misc }

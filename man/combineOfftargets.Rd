\name{combineOfftargets}
\alias{combineOfftargets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine Offtargets
}
\description{
Merge offtargets from different samples
}
\usage{
combineOfftargets(offtarget.folder, sample.name, 
    offtarget.filename = "offTargetsInPeakRegions.xls", 
    common.col = c("targetSeqName", "chromosome", 
        "offTargetStrand", "offTarget_Start", 
        "offTarget_End", "gRNAPlusPAM", "offTarget_sequence",
         "n.mismatch", "guideAlignment2OffTarget", 
         "predicted_cleavage_score"), 
    exclude.col = 2, outputFileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{offtarget.folder}{
   offtarget summary output folder created in GUIDEseqAnalysis function
}
  \item{sample.name}{
   Sample names to be used as part of the column names
   in the final output file
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
  columns to be excluded before merging. Default to 2 (second column of 
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
if(interactive())
{
    offtarget.folder 
    sample.name
    outputFileName
    mergedOfftargets <- 
        combineOfftargets(offtarget.folder = offtarget.folder, 
        sample.name = sample.name, 
        outputFileName = outputFileName) 
}
}
\keyword{ misc }

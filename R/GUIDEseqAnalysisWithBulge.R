GUIDEseqAnalysisWithBulge <- function(alignment.inputfile,
                          umi.inputfile,
                          BSgenomeName,
                          gRNA.file,
                          outputDir,
                          keepPeaksInBothStrandsOnly = FALSE,
                          txdb,
                          orgAnn,
                          PAM.size = 3L,
                          gRNA.size = 20L,
                          PAM = "NGG",
                          ...)
{
  res.mismatch <- GUIDEseqAnalysis(alignment.inputfile = alignment.inputfile,
                          umi.inputfile = umi.inputfile,
                          BSgenomeName = BSgenomeName,
                          gRNA.file = gRNA.file,
                          keepPeaksInBothStrandsOnly = keepPeaksInBothStrandsOnly,
                          txdb = txdb,
                          orgAnn = orgAnn,
                          PAM.size = PAM.size,
                          gRNA.size = gRNA.size,
                          overlap.gRNA.positions = overlap.gRNA.positions,
                          outputDir = outputDir,
                          PAM.location = PAM.location,
                          mat,
                          mismatch.activity.file,
                          ...)
  if (class(gRNA.file) != "DNAStringSet")
  {
    gRNA <- readDNAStringSet(gRNA.file, use.names = TRUE)
    gRNA.name <- names(gRNA)
    gRNA <- as.character(gRNA)
  }
  else
  {
    gRNA.name <- names(gRNA.file)
    gRNA <- as.character(gRNA.file)
  }
  gRNA <- substr(gRNA, 1, gRNA.size)

   res.bulge <- offTargetAnalysisWithBulge(gRNA = gRNA, gRNA.name = gRNA.name,
                                 peaks = "gRNA-PlusMinusPeaksMerged.bed",
                                 BSgenomeName = BSgenomeName,
                                 ...)
  offtargets.mismatch <- res.mismatch[[1]]
  write.table(offtargets.mismatch,
                file = file.path(outputDir,"offTargetsInPeakRegions.xls"),
                sep="\t", row.names = FALSE)

  write.table(res.bulge,
             file = file.path(outputDir,"offTargetsWithBulgeInPeakRegions.xls"),
             sep="\t", row.names = FALSE)

  list(offtargets.mismatch.only = res.mismatch,
       offtargets.with.bulge = res.bulge)
}

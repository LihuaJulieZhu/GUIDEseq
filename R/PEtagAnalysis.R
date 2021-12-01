PEtagAnalysis <- function(alignment.inputfile,
    umi.inputfile, 
    BSgenomeName,
    gRNA.file,
    outputDir,
    keepPeaksInBothStrandsOnly = FALSE,
    txdb, 
    orgAnn, 
    PAM.size = 3L,
    gRNA.size = 20L,
    overlap.gRNA.positions = c(17,18),
    PAM.location = "3prime",
    PBS.len = 10L,
    HA.len = 7L,
    ...) 
{
   res <- GUIDEseqAnalysis(alignment.inputfile = alignment.inputfile,
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
               ...)
   if (PAM.location == "3prime")
   {
      offtargets <- res[[1]]
      baseBeforegRNA <- PBS.len - overlap.gRNA.positions[1]
      baseAfterPAM <- HA.len - PAM.size - gRNA.size + overlap.gRNA.positions[1]

      if (dim(offtargets)[1] > 0)
      {
        chr <- as.character(offtargets$chromosome)
        strand <- as.character(offtargets$offTargetStrand)
        Start <- ifelse(strand=="-",
              as.numeric(as.character( offtargets$offTarget_Start)) - baseAfterPAM,
              as.numeric(as.character( offtargets$offTarget_Start)) - baseBeforegRNA)
        End <- ifelse(strand=="-",
              as.numeric(as.character( offtargets$offTarget_End)) + as.numeric(baseBeforegRNA),
              as.numeric(as.character( offtargets$offTarget_End)) + as.numeric(baseAfterPAM))
       }
       starts <- unlist(apply(cbind(Start,1), 1, max))
       ends <- unlist(apply(cbind(End, seqlengths(BSgenomeName)[chr]), 1,min))
       extendedSequence <- getSeq(BSgenomeName, names = chr, start = starts,
             end = ends, strand = strand, width = NA, as.character = TRUE)
       PBS <- substr(extendedSequence, 1,  PBS.len)
       HA <- substr(extendedSequence, PBS.len + 1, PBS.len + HA.len)   
       offtargets <- cbind(offtargets, PBS = PBS, HAseq = HA)
       res[[1]] <- offtargets
       write.table(offtargets,
          file = file.path(outputDir,"offTargetsInPeakRegions.xls"),
          sep="\t", row.names = FALSE)      
  }
  res  
}

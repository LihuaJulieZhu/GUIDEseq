test_that("GUIDEseqAnalysis", {
     library(TxDb.Hsapiens.UCSC.hg19.knownGene)
     library(org.Hs.eg.db)
     if("BSgenome.Hsapiens.UCSC.hg38" %in% (.packages()))
     {
          detach("package:BSgenome.Hsapiens.UCSC.hg38", unload=TRUE)
     }
     library(BSgenome.Hsapiens.UCSC.hg19)
     umiFile <- system.file("extdata", "UMI-HEK293_site4_chr13.txt",
           package = "GUIDEseq")
     alignFile <- system.file("extdata","bowtie2.HEK293_site4_chr13.sort.bam" ,
            package = "GUIDEseq")
     gRNA.file <- system.file("extdata","gRNA.fa", package = "GUIDEseq")
    cat("GUIDEseqAnalysis without bulge. \n")
    expect_warning(
      GS.res <- GUIDEseqAnalysis(
        alignment.inputfile = alignFile,
        umi.inputfile = umiFile,
        gRNA.file = gRNA.file,
        orderOfftargetsBy = "peak_score",
        descending = TRUE,
        keepTopOfftargetsBy = "predicted_cleavage_score",
        scoring.method = "CFDscore",
        BSgenomeName = Hsapiens,
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
        orgAnn = org.Hs.egSYMBOL,
        outputDir = "GUIDEseqTestResults",
        min.reads = 80, n.cores.max = 1,
        keepPeaksInBothStrandsOnly = FALSE)
    )

    GS.res <- GUIDEseqAnalysis(
        alignment.inputfile = alignFile,
        umi.inputfile = umiFile,
        gRNA.file = gRNA.file,
        orderOfftargetsBy = "peak_score",
        descending = TRUE,
        keepTopOfftargetsBy = "predicted_cleavage_score",
        scoring.method = "CFDscore",
        BSgenomeName = Hsapiens,
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
        orgAnn = org.Hs.egSYMBOL,
        outputDir = "GUIDEseqTestResults",
        min.reads = 80, n.cores.max = 1,
        keepPeaksInBothStrandsOnly = FALSE,
        resume = TRUE)
    
    cat("GUIDEseqAnalysis with offtargets containing no bulge and with bulge. \n")
    expect_warning(GSBulge.res1 <- GUIDEseqAnalysis(
        alignment.inputfile = alignFile,
        umi.inputfile = umiFile,
        gRNA.file = gRNA.file,
        orderOfftargetsBy = "peak_score",
        descending = TRUE,
        keepTopOfftargetsBy = "predicted_cleavage_score",
        scoring.method = "CFDscore",
        BSgenomeName = Hsapiens,
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
        orgAnn = org.Hs.egSYMBOL,
        outputDir = "GUIDEseqTestResults2",
        min.reads = 80, n.cores.max = 1,
        keepPeaksInBothStrandsOnly = FALSE,
        includeBulge = TRUE
    ))

     GSBulge.res1 <- GUIDEseqAnalysis(
        alignment.inputfile = alignFile,
        umi.inputfile = umiFile,
        gRNA.file = gRNA.file,
        orderOfftargetsBy = "peak_score",
        descending = TRUE,
        keepTopOfftargetsBy = "predicted_cleavage_score",
        scoring.method = "CFDscore",
        BSgenomeName = Hsapiens,
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
        orgAnn = org.Hs.egSYMBOL,
        outputDir = "GUIDEseqTestResults2",
        min.reads = 80, n.cores.max = 1,
        keepPeaksInBothStrandsOnly = FALSE,
        includeBulge = TRUE, resume = TRUE
    )
    res2 <- GSBulge.res1$offTargets
    res1 <- GS.res$offTargets
    OT <- "chr13:+:39262912:39262934"
    expect_equal(res1[res1$offTarget == OT,]$guideAlignment2OffTarget,
      res2[res2$offTarget == OT,]$guideAlignment2OffTarget)
   
    expect_equal(res1[res1$offTarget == OT,]$offTargets$predicted_cleavage_score,
      res2[res2$offTarget == OT,]$offTargets$predicted_cleavage_score)

    OT <- "chr13:-:27629404:27629426"
    expect_equal(res1[res1$offTarget == OT,]$guideAlignment2OffTarget,
      res2[res2$offTarget == OT,]$guideAlignment2OffTarget)
   
    expect_equal(res1[res1$offTarget == OT,]$offTargets$predicted_cleavage_score,
      res2[res2$offTarget == OT,]$offTargets$predicted_cleavage_score)
    
    OT <- "chr13+:30964712:30964730" 
    expect_equal(res2[res2$offTarget == OT,]$guideAlignment2OffTarget,
       "A.....T-.....TTTT.TT")

   OT <- "chr13+:67161908:67161923"
    expect_equal(res2[res2$offTarget == OT,]$guideAlignment2OffTarget,
       "A...A.A.AA.A.A..-..A")

   cat("GUIDEseqAnalysis with bulge containing offtargets where no offtargets contain no bulge. \n")


   expect_error(GSBulge.res <- suppressWarnings(GUIDEseqAnalysis(
      alignment.inputfile = alignFile,
      umi.inputfile = umiFile,
      gRNA.file = DNAStringSet("GGCACTGCGGCGGAGGTGGA"),
      orderOfftargetsBy = "peak_score",
      descending = TRUE,
      keepTopOfftargetsBy = "predicted_cleavage_score",
      scoring.method = "CFDscore",
      BSgenomeName = Hsapiens,
      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
      orgAnn = org.Hs.egSYMBOL,
      outputDir = "GUIDEseqTestResults3",
      min.reads = 60, n.cores.max = 1,
      keepPeaksInBothStrandsOnly = FALSE,
      includeBulge = TRUE)))

   if(exists("GSBulge.res"))
   {
      print(GSBulge.res$offTargets)

     expect_equal(GSBulge.res$offTargets$guideAlignment2OffTarget,
        c("A...G......T^.......T", ".......G..T^........G"))
     expect_equal(as.numeric(GSBulge.res$offTargets$predicted_cleavage_score),
        c(0.000213, 0.000259))
#    expect_equal(GSBulge.res$offTargets$mismatch.type,
#             c("rG:dT,rC:dC,rG:dA,rA:dA", "rC:dC,rC:dA,rA:dC"))
     expect_equal(as.character(GSBulge.res$offTargets$DNA.bulge),
             c("A", "U"))
#    expect_equal(GSBulge.res$offTargets$pos.mismatch,
#             list(c(1,5,12,20), c(8, 11, 20 )))
     expect_equal(as.character(GSBulge.res$offTargets$mismatch.distance2PAM),
             c("20,16,9,1", "13,10,1"))
     expect_equal(as.numeric(as.character(GSBulge.res$offTargets$pos.DNA.bulge)), c(13,12))
  }
})

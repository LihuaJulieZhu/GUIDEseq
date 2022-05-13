test_that("GUIDEseqAnalysis", {
     library(TxDb.Hsapiens.UCSC.hg19.knownGene)
     library(org.Hs.eg.db)
     library(BSgenome.Hsapiens.UCSC.hg19)
     umiFile <- system.file("extdata", "UMI-HEK293_site4_chr13.txt",
           package = "GUIDEseq")
     alignFile <- system.file("extdata","bowtie2.HEK293_site4_chr13.sort.bam" ,
            package = "GUIDEseq")
     gRNA.file <- system.file("extdata","gRNA.fa", package = "GUIDEseq")
     expect_warning(GS.res <- GUIDEseqAnalysis(
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
        outputDir = "PEtagTestResults",
        min.reads = 80, n.cores.max = 1,
        keepPeaksInBothStrandsOnly = FALSE
    ))
     GSBulge.res <- GUIDEseqAnalysis(
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
        outputDir = "PEtagTestResults",
        min.reads = 80, n.cores.max = 1,
        keepPeaksInBothStrandsOnly = FALSE,
        includeBulge = TRUE
    )   

    expect_equal(GSBulge.res$offTargets$guideAlignment2OffTarget,
      GS.res$offTargets$guideAlignment2OffTarget)
  
    expect_equal(GSBulge.res$offTargets$predicted_cleavage_score, 
      GS.res$offTargets$predicted_cleavage_score)
})

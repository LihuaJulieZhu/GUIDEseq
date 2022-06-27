test_that("PEtagAnalysis", {
     library(TxDb.Hsapiens.UCSC.hg19.knownGene)
     library(org.Hs.eg.db)
     if("BSgenome.Hsapiens.UCSC.hg38" %in% (.packages()))
     {
          detach("package:BSgenome.Hsapiens.UCSC.hg38", unload=TRUE)
     }
     library(BSgenome.Hsapiens.UCSC.hg19)
     PET.res.truth <- readRDS(system.file("extdata",
            "PEtagAnalysis_testthat.RDS", package = "GUIDEseq"))
     umiFile <- system.file("extdata", "UMI-HEK293_site4_chr13.txt",
           package = "GUIDEseq")
     alignFile <- system.file("extdata","bowtie2.HEK293_site4_chr13.sort.bam" ,
            package = "GUIDEseq")
     gRNA.file <- system.file("extdata","gRNA.fa", package = "GUIDEseq")
     PET.res <- PEtagAnalysis(
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
        PBS.len = 10L,
        HA.len = 7L, min.read.coverage = 2L
    )
  expect_equal(sum(PET.res$offTargets$peak_score),
               sum(PET.res.truth$offTargets$peak_score))
  expect_equal(max(PET.res$offTargets$peak_score),
               max(PET.res.truth$offTargets$peak_score))
  expect_equal(as.vector(PET.res$offTargets[PET.res$offTargets$offTarget ==
                                    "chr13:+:39262912:39262934", 1:20]),
               as.vector(PET.res.truth$offTargets[PET.res.truth$offTargets$offTarget ==
                                          "chr13:+:39262912:39262934", 1:20]))
})

test_that("PEtagAnalysis", {
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
            HA.len = 7L
            )

  expect_equal(PET.res, PET.res.truth)
})

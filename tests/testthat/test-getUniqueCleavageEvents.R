test_that("getUniqueCleavageEvents", {
  cleavages.truth <- readRDS(system.file("extdata", "cleavages_testthat.RDS",
                                         package = "GUIDEseq"))
  umiFile <- system.file("extdata", "UMI-HEK293_site4_chr13.txt",
                         package = "GUIDEseq")
  alignFile <- system.file("extdata","bowtie2.HEK293_site4_chr13.sort.bam" ,
                           package = "GUIDEseq")
  cleavages <- getUniqueCleavageEvents(
    alignment.inputfile = alignFile , umi.inputfile = umiFile,
    n.cores.max = 1)
  expect_equal(cleavages, cleavages.truth)
})

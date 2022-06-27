test_that("getUniqueCleavageEvents", {
  cleavages.truth <- readRDS(system.file("extdata", "cleavages_testthat.RDS",
                                         package = "GUIDEseq"))
  umiFile <- system.file("extdata", "UMI-HEK293_site4_chr13.txt",
                         package = "GUIDEseq")
  alignFile <- system.file("extdata","bowtie2.HEK293_site4_chr13.sort.bam" ,
                           package = "GUIDEseq")
  cleavages <- getUniqueCleavageEvents(
   alignment.inputfile = alignFile , umi.inputfile = umiFile,
   n.cores.max = 1,
   min.umi.count = 5L, max.umi.count = 100000L, min.read.coverage = 1L)
  expect_equal(length(cleavages$cleavage.gr),
               length(cleavages.truth$cleavage.gr))
  expect_equal(sum(cleavages$umi.count.summary$n),
               sum(cleavages.truth$umi.count.summary$n))
  expect_equal(max(cleavages$umi.count.summary$n),
               max(cleavages.truth$umi.count.summary$n))
  expect_equal(median(cleavages$umi.count.summary$n),
               median(cleavages.truth$umi.count.summary$n))
  expect_equal(sd(cleavages$umi.count.summary$n),
               sd(cleavages.truth$umi.count.summary$n))
})

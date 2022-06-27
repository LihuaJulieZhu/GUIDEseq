test_that(".getPAMInd works", {
  #sequences <- readRDS( "~/DropboxUmass/Bioconductor/Trunk/GUIDEseq/inst/extdata/1450peakSeq.RDS")

  #sequences <- readRDS(system.file("extdata", "1450peakSeq.RDS",
#        package = "GUIDEseq"))
  #gRNA <- "TGCTTGGTCGGCACTGATAG"

  sequences <- c("TAATGTTGTTTTCGTTTCTAGGAAGG",
       "TCCTTCCATATATATATATATATATATATATAT")
  PAM.pattern <- "NGG"
  rev.pattern <- "CCN"
  ind.f <- .getPAMInd(PAM.pattern = PAM.pattern, sequences = sequences)
  ind.r <- .getPAMInd(PAM.pattern = rev.pattern , sequences = sequences)
  testthat::expect_equal(as.numeric(ind.f[[1]]), c(20,24))
  testthat::expect_equal(as.numeric(ind.f[[2]]), -1)
  testthat::expect_equal(as.numeric(ind.r[[1]]), -1)
  testthat::expect_equal(as.numeric(ind.r[[2]]), c(2,6))
})

test_that(".getMatchedPAMInd works", {
  sequences <- c("TAATGTTGTTTTCGTTTCTAGGAAGG",
       "TCCTTCCATATATATATATATATATATATATAT")
  PAM.pattern <- "NGG"
  rev.pattern <- "CCN"
  ind.f <- .getPAMInd(PAM.pattern = PAM.pattern, sequences = sequences)
  ind.r <- .getPAMInd(PAM.pattern = rev.pattern , sequences = sequences)

  PAM.location = "3prime"
  PAM.size = 3L
  gRNA.size = 20L
  expect_warning(temp <- .getMatchedPAMInd(inds=  ind.f[[1]],
                         direction = "forward",
                           sequences = sequences[1],
                           PAM.size = PAM.size,
                           gRNA.size = gRNA.size,
                           PAM.location = PAM.location))
  ind.start <- temp[[1]]
  ind.end <- temp[[2]]
  testthat::expect_equal(as.numeric(ind.start),c(1,4))
  testthat::expect_equal(as.numeric(ind.end),c(19,23))

  temp <- .getMatchedPAMInd(inds=  ind.f[[2]],
                         direction = "forward",
                           sequences = sequences[2],
                           PAM.size = PAM.size,
                           gRNA.size = gRNA.size,
                           PAM.location = PAM.location)
  ind.start <- temp[[1]]
  ind.end <- temp[[2]]
  testthat::expect_equal(as.numeric(ind.start),1)
  testthat::expect_equal(as.numeric(ind.end),1)

  temp <- .getMatchedPAMInd(inds=  ind.r[[1]],
                         direction = "reverse",
                           sequences = sequences[1],
                           PAM.size = PAM.size,
                           gRNA.size = gRNA.size,
                           PAM.location = PAM.location)

  ind.start2 <- temp[[1]]
  ind.end2 <- temp[[2]]
  testthat::expect_equal(ind.start2,1)
  testthat::expect_equal(ind.end2,1)

  temp <- .getMatchedPAMInd(inds=  ind.r[[2]],
                         direction = "reverse",
                           sequences = sequences[2],
                           PAM.size = PAM.size,
                           gRNA.size = gRNA.size,
                           PAM.location = PAM.location)

  ind.start2 <- temp[[1]]
  ind.end2 <- temp[[2]]
  testthat::expect_equal(ind.start2, c(5,9))
  testthat::expect_equal(ind.end2, c(24, 28))
})

test_that(".getSubSeq works", {
  sequences <- c("TAATGTTGTTTTCGTTTCTAGGAAGG",
       "TCCTTCCATATATATATATATATATATATATAT")
  PAM.pattern <- "NGG"
  rev.pattern <- "CCN"
  PAM.size <- 3L
  gRNA.size <- 20L
  PAM.location <- "3prime"
  ind.f <- .getPAMInd(PAM.pattern = PAM.pattern, sequences = sequences)
  ind.r <- .getPAMInd(PAM.pattern = rev.pattern , sequences = sequences)

  expect_warning(seq.f <- do.call(rbind, lapply(1:length(sequences),
     function(i) {
       .getSubSeq(ind.f[[i]],
                  sequences[i],
                  direction = "forward", PAM.location = PAM.location,
                  gRNA.size = gRNA.size, PAM.size = PAM.size,
                  max.DNA.bulge = 0L)
   })))

  testthat::expect_equal(seq.f[,5], c("TAATGTTGTTTTCGTTTCT",
          "TGTTGTTTTCGTTTCTAGGA"))
  testthat::expect_equal(substr(sequences[1], 4, 26),
      as.character(seq.f[2,7]))
  testthat::expect_equal(substr(sequences[1], 1, 22),
      as.character(seq.f[1,7]))

  testthat::expect_equal(seq.f[,6], c("AGG",
                                    "AGG"))
  testthat::expect_equal(seq.f[,3], c("1", "4"))
  testthat::expect_equal(seq.f[,4], c("22", "26"))
  testthat::expect_equal(seq.f[,2], c("+", "+"))
  testthat::expect_equal(seq.f[,1], rep("TAATGTTGTTTTCGTTTCTAGGAAGG", 2))


  seq.r <-  do.call(rbind, lapply(1:length(sequences),
      function(i) {
        .getSubSeq(ind.r[[i]],
                   sequences[i], direction = "reverse",
                   PAM.location = PAM.location,
                   gRNA.size = gRNA.size, PAM.size = PAM.size,
                  max.DNA.bulge = 0L)
   }))

   testthat::expect_equal(seq.r[,1],
       rep("TCCTTCCATATATATATATATATATATATATAT",2))
   testthat::expect_equal(seq.r[,2], c("-", "-"))

  testthat::expect_equal(as.character(reverseComplement(
    DNAString(substr(sequences[2], 2, 24)))),
    as.character(seq.r[1,7]))

  testthat::expect_equal(as.character(reverseComplement(
    DNAString(substr(sequences[2], 6, 28)))),
    as.character(seq.r[2,7]))

  testthat::expect_equal(seq.r[,6], c("AGG","TGG"))

  testthat::expect_equal(seq.r[,3], c("2", "6"))
  testthat::expect_equal(seq.r[,4], c("24", "28"))
})


test_that("getMaskedSeq works with max.DNA.bulge of 2", {
  sequences <- c("TAATGTTGTTTTCGTTTCTAGGAAGG",
       "TCCTTCCATATATATATATATATATATATATAT")
  PAM.location <- "3prime"
  PAM.pattern <- "NGG$"
  expect_warning(sub.seq <- .PAMpatternSearch(PAM.pattern, sequences,
                              PAM.location = PAM.location,
                              PAM.size = 3,
                              gRNA.size = 20,
                              max.DNA.bulge = 2L))

   testthat::expect_equal(sub.seq[,1],
       c(rep(sequences[1], 2),
         rep(sequences[2],2)))

   testthat::expect_equal(sub.seq[,2],
       c(rep("+", 2),
         rep("-",2)))

   testthat::expect_equal(as.numeric(sub.seq[,3]), c(1,2,2,6))
   testthat::expect_equal(as.numeric(sub.seq[,4]), c(22,26,26,30))
   testthat::expect_equal(sub.seq[,6], c(rep("AGG", 3), "TGG"))

   testthat::expect_equal(sub.seq[,7], c("TAATGTTGTTTTCGTTTCTAGG",
     "AATGTTGTTTTCGTTTCTAGGAAGG",
     as.character(reverseComplement(DNAString(
       substr("TCCTTCCATATATATATATATATATATATATAT", 2, 26)))),
     as.character(reverseComplement(DNAString(
       substr("TCCTTCCATATATATATATATATATATATATAT", 6, 30))))))
})

test_that(".getSubSeq allowing max.DNA.bulge to 2 works", {
  sequences <- c("TAATGTTGTTTTCGTTTCTAGGAAGG",
                 "TCCTTCCATATATATATATATATATATATATAT")
  PAM.pattern <- "NGG"
  rev.pattern <- "CCN"
  PAM.size <- 3L
  gRNA.size <- 20L
  PAM.location = "3prime"
  ind.f <- .getPAMInd(PAM.pattern = PAM.pattern, sequences = sequences)
  ind.r <- .getPAMInd(PAM.pattern = rev.pattern , sequences = sequences)

  max.DNA.bulge <- 2L

  expect_warning(seq.f <- do.call(rbind, lapply(1:length(sequences),
                                 function(i) {
                                   .getSubSeq(ind.f[[i]],
                                              sequences[i],
                                              direction = "forward",
                                              PAM.location = PAM.location,
                                              gRNA.size = gRNA.size,
                                              PAM.size = PAM.size,
                                              max.DNA.bulge = max.DNA.bulge)
                                 })))

  testthat::expect_equal(seq.f[,5], c("TAATGTTGTTTTCGTTTCT",
                                      "AATGTTGTTTTCGTTTCTAGGA"))
  testthat::expect_equal(substr(sequences[1], 2, 26),
                         as.character(seq.f[2,7]))
  testthat::expect_equal(substr(sequences[1], 1, 22),
                         as.character(seq.f[1,7]))

  testthat::expect_equal(seq.f[,6], c("AGG",
                                      "AGG"))
  testthat::expect_equal(seq.f[,3], c("1", "2"))
  testthat::expect_equal(seq.f[,4], c("22", "26"))
  testthat::expect_equal(seq.f[,2], c("+", "+"))
  testthat::expect_equal(seq.f[,1], rep("TAATGTTGTTTTCGTTTCTAGGAAGG", 2))


  seq.r <-  do.call(rbind, lapply(1:length(sequences),
                                  function(i) {
                                    .getSubSeq(ind.r[[i]],
                                               sequences[i], direction = "reverse",
                                               PAM.location = PAM.location,
                                               gRNA.size = gRNA.size, PAM.size = PAM.size,
                                               max.DNA.bulge = max.DNA.bulge)
                                  }))

  testthat::expect_equal(seq.r[,1],
                         rep("TCCTTCCATATATATATATATATATATATATAT",2))
  testthat::expect_equal(seq.r[,2], c("-", "-"))

  testthat::expect_equal(as.character(reverseComplement(
    DNAString(substr(sequences[2], 2, 26)))),
    as.character(seq.r[1,7]))

  testthat::expect_equal(as.character(reverseComplement(
    DNAString(substr(sequences[2], 6, 30)))),
    as.character(seq.r[2,7]))

  testthat::expect_equal(seq.r[,6], c("AGG","TGG"))

  testthat::expect_equal(seq.r[,3], c("2", "6"))
  testthat::expect_equal(seq.r[,4], c("26", "30"))

  cat("test sequence not long enough for allowing large stretch of bulge")
  max.DNA.bulge <- 8
  seq.r <-  do.call(rbind, lapply(1:length(sequences),
                                  function(i) {
                                    .getSubSeq(ind.r[[i]],
                                               sequences[i], direction = "reverse",
                                               PAM.location = PAM.location,
                                               gRNA.size = gRNA.size, PAM.size = PAM.size,
                                               max.DNA.bulge = max.DNA.bulge)
                                  }))

  testthat::expect_equal(seq.r[,6], c("AGG","TGG"))

  testthat::expect_equal(seq.r[,3], c("2", "6"))
  testthat::expect_equal(seq.r[,4], c("32", "33"))


  testthat::expect_equal(seq.r[,1],
                         rep("TCCTTCCATATATATATATATATATATATATAT",2))
  testthat::expect_equal(seq.r[,2], c("-", "-"))

  expect_warning(seq.f <- do.call(rbind, lapply(1:length(sequences),
                                 function(i) {
                                   .getSubSeq(ind.f[[i]],
                                              sequences[i],
                                              direction = "forward",
                                              PAM.location = PAM.location,
                                              gRNA.size = gRNA.size,
                                              PAM.size = PAM.size,
                                              max.DNA.bulge = max.DNA.bulge)
                                 })))

  testthat::expect_equal(seq.f[,6], c("AGG",
                                      "AGG"))
  testthat::expect_equal(seq.f[,3], c("1", "1"))
  testthat::expect_equal(seq.f[,4], c("22", "26"))
  testthat::expect_equal(seq.f[,2], c("+", "+"))
  testthat::expect_equal(seq.f[,1], rep("TAATGTTGTTTTCGTTTCTAGGAAGG", 2))
  testthat::expect_equal(seq.f[,5], c("TAATGTTGTTTTCGTTTCT",
                                      "TAATGTTGTTTTCGTTTCTAGGA"))
})


test_that(".getSubSeq with default max.DNA.bulge = 2 works for both strands", {
  sequences <- c("TAATGTTGTTTTCGTTTCTAGGAAGG",
                 "TCCTTCCATATATATATATATATATATATATAT")
  PAM.location <- "3prime"
  PAM.pattern <- "NGG$"

  expect_warning(sub.seq <- .PAMpatternSearch(PAM.pattern, sequences,
                               PAM.location = "3prime",
                               PAM.size = 3,
                               gRNA.size = 20))

  testthat::expect_equal(sub.seq[,1],
                         c(rep(sequences[1], 2),
                           rep(sequences[2],2)))

  testthat::expect_equal(sub.seq[,2],
                         c(rep("+", 2),
                           rep("-",2)))

  testthat::expect_equal(as.numeric(sub.seq[,3]), c(1,2,2,6))
  testthat::expect_equal(as.numeric(sub.seq[,4]), c(22,26,26,30))
  testthat::expect_equal(sub.seq[,6], c(rep("AGG", 3), "TGG"))

  testthat::expect_equal(sub.seq[,7], c("TAATGTTGTTTTCGTTTCTAGG",
                                        "AATGTTGTTTTCGTTTCTAGGAAGG",
                                        "TATATATATATATATATATGGAAGG",
                                      "TATATATATATATATATATATATGG"))
})

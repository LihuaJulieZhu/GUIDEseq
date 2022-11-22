# pa.f <- readRDS(file = system.file("extdata", "pa.f.RDS",
#                                    package = "GUIDEseq"))
# pa.r <- readRDS(file = system.file("extdata", "pa.r.RDS",
#                                    package = "GUIDEseq"))
subjects2 <- readRDS(file = system.file("extdata", "subjects2.RDS",
                                         package = "GUIDEseq"))

# hg19
gRNA <- "TTGCTTTTATCACAGGCTCC"

pa.f<- lapply(1:length(subjects2), function(i) {
  pairwiseAlignment(pattern = gRNA,
                    subject = as.character(subjects2[i]),
                    type = "global-local",
                    scoreOnly = FALSE)
})

pa.r <- lapply(1:length(subjects2), function(i) {
  pairwiseAlignment(pattern = gRNA,
                    subject = as.character(reverseComplement(subjects2[i])),
                    type = "global-local",
                    scoreOnly = FALSE)
})

#saveRDS(pa.f, file = "~/Dropbox (UMass Medical School)/Bioconductor/Trunk/GUIDEseq/inst/extdata/pa.f.RDS")
#saveRDS(pa.r, file = "~/Dropbox (UMass Medical School)/Bioconductor/Trunk/GUIDEseq/inst/extdata/pa.r.RDS")


test_that("getBestAlnInfo plus strand bulge in gRNA works", {
  i = 16
  temp <- getBestAlnInfo(subjects2[i], pa.f[[i]], pa.r[[i]])
  expect_equal(
    as.character(reverseComplement(DNAString(substr(as.character(subjects2[i]),
                        temp$offTarget_Start, temp$offTarget_End)))),
    "TTCCATTCAATGATTCCATT")
  expect_equal(temp$offTargetStrand, "-")
  expect_equal(temp$guideAlignment2OffTarget, "..C.A..---..AT.AT...")
  expect_equal(temp$n.mismatch, 6)
  expect_equal(temp$n.insertion, 3)
  expect_equal(temp$n.deletion, 0)
  expect_equal(temp$pos.insertion, c(8,9,10))
  expect_equal(temp$pos.mismatch, c(3, 5, 13, 14, 16, 17))
})

test_that("getBestAlnInfo plus strand bulge in gRNA works", {
  i = 123
  temp <- getBestAlnInfo(subjects2[i], pa.f[[i]], pa.r[[i]])
  expect_equal(
    as.character(substr(as.character(subjects2[i]),
                        temp$offTarget_Start, temp$offTarget_End)),
    "ATGGTTTCAAAACTGCTCCATC")

  expect_equal(temp$offTarget, "chrX+:61721170:61721195chrX-:61721140:61721165")
  expect_equal(as.character(temp$PAM.sequence),"ATC")
  expect_equal(temp$offTarget_sequence, "ATGGTTTCAAAACTG-CTCCATC")
  expect_equal(temp$offTargetStrand, "+")
  expect_equal(temp$offTarget_End, 57)
  expect_equal(temp$offTarget_Start, 36)
  expect_equal(temp$guideAlignment2OffTarget, "A..G...C.AA..T.-....")
  expect_equal(temp$n.PAM.mismatch, 2)
  expect_equal(temp$mismatch.distance2PAM, "20,17,13,11,10,7")
  expect_equal(as.numeric(temp$pos.insertion), 16)
  expect_equal(temp$pos.mismatch, c(1, 4, 8, 10, 11, 14))
  # n.mismatch = mismatches + indels
  expect_equal(temp$n.mismatch, 6)

  expect_equal(temp$n.deletion, 0)
  expect_equal(temp$n.insertion, 1)
})

test_that("getBestAlnInfo minus strand bulge in gRNA works", {
  i = 66
  temp <- getBestAlnInfo(subjects2[i], pa.f[[i]], pa.r[[i]])
  expect_equal(
    as.character(reverseComplement(DNAString(substr(subjects2[i],
                        temp$offTarget_Start, temp$offTarget_End)))),
    "TTGCATTATTTCAGCATCTCAT")

  expect_equal(temp$offTarget,
                         "chr2+:221535820:221535845chr2-:221535730:221535755"
  )

  expect_equal(temp$offTargetStrand, "-")
  expect_equal(as.character(temp$PAM.sequence),"CAT")
  expect_equal(temp$offTarget_sequence, "TTGCATT-ATTTCAGCATCTCAT")

  expect_equal(temp$offTarget_End, 71)
  expect_equal(temp$offTarget_Start, 50)

  expect_equal(temp$guideAlignment2OffTarget, "....A..-..TT...CA..T")
  expect_equal(temp$n.PAM.mismatch, 2)
  expect_equal(temp$mismatch.distance2PAM, "16,10,9,5,4,1") ###
  expect_equal(as.numeric(temp$pos.insertion), 8)
  expect_equal(temp$pos.mismatch, c(5,11,12,16,17,20))
  # n.mismatch = mismatches + indels
  expect_equal(temp$n.mismatch, 6)
  expect_equal(temp$pos.insertion, 8)

  expect_equal(temp$n.insertion, 1)
})

test_that("getBestAlnInfo plus strand without indel works", {
  i = 1
  #pattern:      TTGCTTTTATCACAGGCTCC
  #subject: [84] GTGTGTTTGGAAACTGCTCC
  temp <- getBestAlnInfo(subjects2[i], pa.f[[i]], pa.r[[i]])
  expected_mismatch_pos <- c(1, 4,5, 9, 10, 11, 13, 14, 15)

  expect_equal(temp$pos.mismatch, expected_mismatch_pos)

  expect_equal(temp$offTarget_End, 106)
  expect_equal(temp$offTarget_Start, 84)

  expect_equal(
    as.character(substr(as.character(subjects2[i]),
                        temp$offTarget_Start, temp$offTarget_End)),
    temp$offTarget_sequence)

  expect_equal(temp$offTarget,names(subjects2[i]))

  expect_equal(temp$offTargetStrand, "+")

  expect_equal(as.character(temp$PAM.sequence),"ATC")

  expect_equal(temp$guideAlignment2OffTarget, "G..TG...GGA.ACT.....")
  expect_equal(temp$n.PAM.mismatch, 2)
  expect_equal(temp$mismatch.distance2PAM, "20,17,16,12,11,10,8,7,6")
   # n.mismatch = mismatches + indels
  expect_equal(temp$n.mismatch, 9)

  expect_equal(temp$n.insertion, 0)
  expect_equal(temp$n.deletion, 0)
})

test_that("getBestAlnInfo minus strand with deletion on gRNA works", {
  i = 96
  #pattern:      GGAGCCTGTGATAAAAGCAA
  #subject: [48] GAAGCCTATGCTAGAAATGG
  temp <- getBestAlnInfo(subjects2[i], pa.f[[i]], pa.r[[i]])

  expect_equal(temp$offTarget_End, 69)
  expect_equal(temp$offTarget_Start, 45)

  expect_equal(
    as.character(substr(as.character(subjects2[i]),
                        temp$offTarget_Start, temp$offTarget_End)),
    as.character(reverseComplement(DNAString(temp$offTarget_sequence))))

  expect_equal(temp$offTarget,names(subjects2[i]))

  expect_equal(temp$offTargetStrand, "-")

  expect_equal(as.character(reverseComplement(DNAString(
                  as.character(temp$PAM.sequence)))),
                         as.character(substr(as.character(subjects2[i]),
                             temp$offTarget_Start, temp$offTarget_Start + 2)))


  expect_equal(temp$guideAlignment2OffTarget, "..C.A..^^..G..T.....T.")
  expect_equal(temp$mismatch.distance2PAM, "18,16,11,8,2")
  expect_equal(temp$pos.mismatch, c(3,5, 10, 13, 19))

  expect_equal(temp$n.PAM.mismatch, 1)
  # n.mismatch = mismatches + indels
  expect_equal(temp$n.mismatch, 5)

  expect_equal(temp$n.deletion, 2)
  expect_equal(temp$n.insertion, 0)
})

#Global-Local PairwiseAlignmentsSingleSubject (1 of 1)
#pattern:      TGCTTGGTCGG-CACTGATAG
#subject: [15] GGCTTGGCCGGGCACTGATTG
test_that("getBestAlnInfo plus strand with bulge in offtargets works", {

   chr14OT <- readRDS(file = system.file(
                             "extdata", "InsertionInOfftarget.RDS",
                                           package = "GUIDEseq"))
  subjects1 <- chr14OT[[3]]
  gRNA <- "TGCTTGGTCGGCACTGATAG"

  pa.f<- lapply(1:length(subjects1), function(i) {
    pairwiseAlignment(pattern = gRNA,
                      subject = as.character(subjects1[i]),
                      type = "global-local",
                      scoreOnly = FALSE)
  })

  pa.r <- lapply(1:length(subjects1), function(i) {
    pairwiseAlignment(pattern = gRNA,
                      subject = as.character(reverseComplement(subjects1[i])),
                      type = "global-local",
                      scoreOnly = FALSE)
  })

  i <- 1

  temp <- getBestAlnInfo(subjects1[i], pa.f[[i]], pa.r[[i]])

  expect_equal(as.numeric(temp$pos.deletion), 10)
  expect_equal(temp$pos.mismatch, c(1,8,19))

  expect_equal(temp$offTarget_End, 38)
  expect_equal(temp$offTarget_Start, 15)

  expect_equal(
    as.character(substr(as.character(subjects1[i]),
                        temp$offTarget_Start, temp$offTarget_End)),
    temp$offTarget_sequence)

  expect_equal(temp$offTarget,names(subjects1[i]))

  expect_equal(temp$offTargetStrand, "+")

  PAM.size <- 3
  expect_equal(
    as.character(temp$PAM.sequence),
    as.character(substr(as.character(subjects1[i]),
                        temp$offTarget_End - PAM.size  + 1, temp$offTarget_End)))

  expect_equal(temp$pos.deletion, 10)

  expect_equal(temp$pos.mismatch, c(1,8,19))
  expect_equal(temp$guideAlignment2OffTarget, "G......C.^.........T.")
  expect_equal(temp$mismatch.distance2PAM, "20,13,2")

  expect_equal(temp$n.PAM.mismatch, 0)
  # n.mismatch = mismatches + indels
  expect_equal(temp$n.mismatch, 3)

  expect_equal(temp$n.insertion, 0)
  expect_equal(temp$n.deletion, 1)
})

test_that("getBestAlnInfo minus strand with perfect match in offtargets works", {

   chr14OT <- readRDS(file = system.file(
     "extdata", "InsertionInOfftarget.RDS",
     package = "GUIDEseq"))

  subjects1 <- chr14OT[[3]]
  gRNA <- "TGCTTGGTCGGCACTGATAG"

  pa.f<- lapply(1:length(subjects1), function(i) {
    pairwiseAlignment(pattern = gRNA,
                      subject = as.character(subjects1[i]),
                      type = "global-local",
                      scoreOnly = FALSE)
  })

  pa.r <- lapply(1:length(subjects1), function(i) {
    pairwiseAlignment(pattern = gRNA,
                      subject = as.character(reverseComplement(subjects1[i])),
                      type = "global-local",
                      scoreOnly = FALSE)
  })

  #saveRDS(list(pa.f = pa.f, pa.r = pa.r, subjects = subjects1), file =
  #         "~/Dropbox (UMass Medical School)/Bioconductor/Trunk/GUIDEseq/inst/extdata/InsertionInOfftarget.RDS")


  #pa.f <- chr14OT[[1]]
  #pa.r <- chr14OT[[2]]

  i <- 2

  temp <- getBestAlnInfo(subjects1[i], pa.f[[i]], pa.r[[i]])

  expect_equal(temp$offTarget_End, 47)
  expect_equal(temp$offTarget_Start, 25)

  expect_equal(
    as.character(substr(as.character(subjects1[i]),
                        temp$offTarget_Start, temp$offTarget_End)),
    as.character(reverseComplement(DNAString(temp$offTarget_sequence))))

  expect_equal(temp$offTarget,names(subjects1[i]))

  expect_equal(temp$offTargetStrand, "-")

  expect_equal(as.character(reverseComplement(DNAString(
    as.character(temp$PAM.sequence)))),
    as.character(substr(as.character(subjects1[i]),
                        temp$offTarget_Start, temp$offTarget_Start + 2)))


  expect_equal(temp$guideAlignment2OffTarget, paste0(rep(".", 20),
                                                               collapse = ""))
  expect_equal(temp$mismatch.distance2PAM, "")

  expect_equal(temp$n.PAM.mismatch, 0)
  # n.mismatch = mismatches + indels
  expect_equal(temp$n.mismatch, 0)

  expect_equal(temp$n.insertion, 0)
  expect_equal(temp$n.deletion, 0)
})

test_that("getBestAlnInfo plus strand with bulge in both guide (first) and OT works", {

   gRNA <- "CTTGCCCCACAGGGCAGTAA"
   peak <- "CTAGCTCTCACAGGAACTAA"
   pa.f<- pairwiseAlignment(pattern = gRNA,
                    subject = peak,
                    type = "global-local",
                    scoreOnly = FALSE)

   pa.r <- pairwiseAlignment(pattern = gRNA,
               subject = as.character(reverseComplement(DNAString(peak))), 
               type = "global-local",
   	       scoreOnly = FALSE)

   temp <- getBestAlnInfo(peak, pa.f, pa.r)
   expect_equal(temp$offTarget_sequence, "CTAGCTCTCACAGGA-ACTAA")
   expect_equal(temp$guideAlignment2OffTarget , "..A..T.^......A-.C...")
   expect_equal(temp$offTargetStrand, "+")
   expect_equal(temp$pos.mismatch, c(3, 6, 14, 17)) 
   expect_equal(temp$n.insertion, 1)
   expect_equal(temp$n.deletion, 1)
   expect_equal(temp$pos.insertion, 15)
   expect_equal(temp$pos.deletion, 8)
})

test_that("getBestAlnInfo plus strand with bulge in both guide (second) and OT (first) works", {
   peak <- "CTTGCCCCACAGGGCAGTAA"
   gRNA <- "CTAGCTCTCACAGGAACTAA"
   pa.f<- pairwiseAlignment(pattern = gRNA,
                    subject = peak,
                    type = "global-local",
                    scoreOnly = FALSE)

   pa.r <- pairwiseAlignment(pattern = gRNA,
               subject = as.character(reverseComplement(DNAString(peak))),
               type = "global-local",
               scoreOnly = FALSE)

   temp <- getBestAlnInfo(peak, pa.f, pa.r)
   expect_equal(temp$offTarget_sequence, "CTTGCCC-CACAGGGCAGTAA")
   expect_equal(temp$offTargetStrand, "+")
   expect_equal(temp$pos.mismatch, c(3, 6, 15, 17))
   expect_equal(temp$guideAlignment2OffTarget, "..T..C.-......G^.G...")
   expect_equal(temp$n.insertion, 1)
   expect_equal(temp$n.deletion, 1)
   expect_equal(temp$pos.insertion, 8)
   expect_equal(temp$pos.deletion, 16)
})

test_that("getBestAlnInfo plus strand with bulge in both guide (second) and OT (first) works", {
   gRNA <- "CTTGCCCCACAGGGCAGTAA"
   peak <- "ATTGCCCAGGCTTAGTGCAGTAGCGC"
   pa.f<- pairwiseAlignment(pattern = gRNA,
                    subject = peak,
                    type = "global-local",
                    scoreOnly = FALSE)

   pa.r <- pairwiseAlignment(pattern = gRNA,
               subject = as.character(reverseComplement(DNAString(peak))),
               type = "global-local",
               scoreOnly = FALSE)

   temp <- getBestAlnInfo(peak, pa.f, pa.r)
   expect_equal(temp$offTarget_sequence, peak)
   expect_equal(temp$offTargetStrand, "+")
   expect_equal(temp$pos.mismatch, c(1, 9, 10, 13, 20))
   expect_equal( temp$guideAlignment2OffTarget, "A......^^^.TT..T......G")
   expect_equal(temp$n.insertion, 0)
   expect_equal(temp$n.deletion, 3)
   expect_equal(temp$pos.deletion, c(8,9,10))
})


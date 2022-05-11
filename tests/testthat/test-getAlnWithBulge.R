# need to remove the source later
source("~/Dropbox (UMass Medical School)/Bioconductor/Trunk/GUIDEseq/R/getAlnWithBulge.R")
source("~/Dropbox (UMass Medical School)/Bioconductor/Trunk/GUIDEseq/R/getBestAlnInfo.R")
source("~/Dropbox (UMass Medical School)/Bioconductor/Trunk/GUIDEseq/R/bulge-internal.R")

detach("package:BSgenome.Hsapiens.UCSC.hg19", unload = TRUE)
library(BSgenome.Hsapiens.UCSC.hg38)
test_that("getAlnWithBulge without restricting PAM.pattern works", {

  #peaks.f <- "~/DropboxUmass/Bioconductor/Trunk/GUIDEseq/inst/extdata/1450-PlusMinusPeaksMerged.bed"
  peaks.f <- "~/DropboxUmass/Bioconductor/Trunk/GUIDEseq/inst/extdata/1450-chr14-chr2-bulge-test.bed"
  gRNA <- "TGCTTGGTCGGCACTGATAG"

  temp <- getAlnWithBulge(gRNA, gRNA.name = "gRNA1450", peaks.withHeader = FALSE,
                          peaks = peaks.f, BSgenomeName = Hsapiens,
  )

  aln.all <- temp$aln.all
  expect_equal(as.character(aln.all$offTarget_sequence),
               c("GGCTTGGCCGGGCACTGATTGCGG", "TGCTTGGTCGGCACTGATAGGGG"))

  expect_equal(aln.all$pos.mismatch[[1]], c(1, 8, 19))
  expect_equal(as.numeric(aln.all$pos.deletion)[1],  12)

  expect_equal(aln.all$guideAlignment2OffTarget[[1]],
               "G......C...^.......T.")

  expect_equal(aln.all$guideAlignment2OffTarget[[2]],
               "....................")

  #featureVectors <- buildFeatureVectorForScoringBulge(aln.all)
})

detach("package:BSgenome.Hsapiens.UCSC.hg38", unload = TRUE)
library("BSgenome.Hsapiens.UCSC.hg19")
test_that("getAlnWithBulge bulge on gRNA offtarget on minus strand without mismatches works", {
  peaks <- DNAStringSet(c("CCTGAGGCTGGGGTGGAGGGGGTC"))
  names(peaks) <- "testMinusBulgeOff"
  gRNA <- substr(as.character(readDNAStringSet(system.file("extdata", "T2.fa",
                                                           package = "CRISPRseek"))),
                 1, 20)
  temp1 <- getAlnWithBulge(gRNA, gRNA.name = "T2",
                           peaks = peaks, BSgenomeName = Hsapiens)

  expect_equal(as.character(temp1$aln.all$guideAlignment2OffTarget),
               "...............^.....")

  expect_true(as.character(temp1$aln.all$offTarget_sequence) ==
                as.character(reverseComplement(peaks)))

  expect_equal(as.numeric(temp1$aln.all$pos.deletion), 16)
})

test_that(" bulge on offtarget with mismacth and on minus strand works", {
  # PAM followed by 2t(19A), 19c(2G), 20t(1A), 6 insertion
  peaks <- DNAStringSet(c("CCTGTGGCTGGGGTGGAGGGGGCT"))
  gRNA <- substr(as.character(readDNAStringSet(system.file("extdata", "T2.fa",
                                                           package = "CRISPRseek"))), 1, 20)

  names(peaks) <- "testMinusBulgeOff"
  temp1 <- getAlnWithBulge(gRNA, gRNA.name = "T2",
                           peaks = peaks, BSgenomeName = Hsapiens)

  expect_equal(as.character(temp1$aln.all$guideAlignment2OffTarget),
               "AG.............^...A.")

  expect_true(as.character(reverseComplement(DNAStringSet(substr(peaks,
                                                                 temp1$aln.all$offTarget_Start, temp1$aln.all$offTarget_End)))) ==
                temp1$aln.all$offTarget_sequence)

  expect_true(as.character(temp1$aln.all$offTarget_sequence) ==
                as.character(reverseComplement(peaks)))

  expect_equal(as.numeric(temp1$aln.all$pos.deletion), 16)

  expect_equal(as.numeric(temp1$aln.all$offTarget_End), 24)
  expect_equal(as.numeric(temp1$aln.all$offTarget_Start), 1)

  expect_equal(as.character(temp1$aln.all$PAM.sequence), "AGG")
  expect_equal(
    as.numeric(unlist(temp1$aln.all$pos.mismatch)),c(1,2,19))
  expect_equal(as.character(temp1$aln.all$offTargetStrand), "-")
})

test_that("bulge on gRNA offtarget on plus strand works", {
  peaks.f <- system.file("extdata", "T2plus100OffTargets.bed",
                         package = "GUIDEseq")
  gRNA <- substr(as.character(readDNAStringSet(system.file("extdata", "T2.fa",
                                                           package = "CRISPRseek"))),
                 1, 20)
  names(gRNA) <- "T2"
  temp <- getAlnWithBulge(gRNA, gRNA.name = "T2",
                          peaks = peaks.f, BSgenomeName = Hsapiens)

  bed <- read.table(system.file("extdata", "T2plus100OffTargets.bed",
                                package = "GUIDEseq"),
                    sep = "\t",
                    header = TRUE)
  merged.bed <- merge(bed, temp$aln.all, by.x = "names", by.y ="offTarget")
  # some sequences do not contain the PAM.pattern
  #expect_equal(nrow(merged.bed), nrow(temp))
  expect_equal(merged.bed$peak_score,merged.bed$totalCount)
  expect_equal(merged.bed$chr.x,merged.bed$chr.y)

  off.names <- c("chr1+:121373900:121373925 chr1-:121373840:121373865",
                 "chr2+:221535820:221535845 chr2-:221535730:221535755",
                 "chr10+:42529870:42529895 chr10-:42529770:42529795")
  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names == off.names[1],]$pos.insertion))
    , 7)
  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[1],]$guideAlignment2OffTarget)),
    "..G.GG-....AAT.T..A.")
  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[1],]$n.mismatch)), 8)
  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names == off.names[1],]$PAM.sequence)),
    "TTG")
  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[1],]$pos.mismatch)),
    c(3, 5, 6, 12, 13, 14, 16, 19))

  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[1],]$offTarget_End)),
    121373887)

  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[1],]$offTarget_Start)),
    121373866)

  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[1],]$chromosome.y)),
    "chr1")

  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[1],]$mismatch.distance2PAM)),
    "18,16,15,9,8,7,5,2")

  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[1],]$n.mismatch)),
    8)

  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[1],]$offTarget_sequence)),
    "GAGCGG-TCCAAATCTCCACTTG")

  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[1],]$offTargetStrand)),
    "+")

  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[3],]$offTarget_Start)),
    42529795)

  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[3],]$offTarget_End)),
    42529816)

  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[3],]$chromosome.y)),
    "chr10")
})

test_that("bulge on gRNA and offtarget on the minus strand  works", {

  cat("Start testing offtarget on the minus strand ...")

  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[2],]$offTarget_sequence)),
    "CACCCCC-CCACCCCACCCCAGG")
  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[2],]$offTargetStrand)),
    "-")

  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[2],]$n.mismatch)),
    3)

  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[2],]$mismatch.distance2PAM)),
    "20,5,2")


  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names == off.names[2],]$pos.insertion))
    , 8)

  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[2],]$guideAlignment2OffTarget)),
    "C......-.......A..C.")

  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[2],]$n.mismatch)), 3)
  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names == off.names[2],]$PAM.sequence)),
    "AGG")
  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[2],]$pos.mismatch)),
    c(1,16,19))



  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[2],]$offTarget_Start)),
    221535826)

  expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[2],]$offTarget_End)),
    221535847)

  expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[2],]$chromosome.y)),
    "chr2")
})

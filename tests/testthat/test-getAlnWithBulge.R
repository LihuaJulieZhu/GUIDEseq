require(CRISPRseek)
require(GUIDEseq)

source("~/Dropbox (UMass Medical School)/Bioconductor/Trunk/GUIDEseq/R/getAlnWithBulge.R")
source("~/Dropbox (UMass Medical School)/Bioconductor/Trunk/GUIDEseq/R/getBestAlnInfo.R")
source("~/Dropbox (UMass Medical School)/Bioconductor/Trunk/GUIDEseq/R/bulge-internal.R")

test_that("getAlnWithBulge without restricting PAM.pattern works", {
  require(BSgenome.Hsapiens.UCSC.hg38)

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
  expect_equal(as.numeric(aln.all$pos.insertion)[1],  12)

  expect_equal(aln.all$guideAlignment2OffTarget[[1]],
               "G......C...^.......T.")

  expect_equal(aln.all$guideAlignment2OffTarget[[2]],
               "....................")

#featureVectors <- buildFeatureVectorForScoringBulge(aln.all)
})

test_that("getAlnWithBulge without restricting PAM.pattern works", {
  ### hg19
  cat("Start testing bulge on gRNA offtarget on minus strand without mismatches...")

  peaks <- DNAStringSet(c("CCTGAGGCTGGGGTGGAGGGGGTC"))
  names(peaks) <- "testMinusBulgeOff"
  gRNA <- substr(as.character(readDNAStringSet(system.file("extdata", "T2.fa",
                                                           package = "CRISPRseek"))),
                 1, 20)
  temp1 <- getAlnWithBulge(gRNA, gRNA.name = "T2",
                           peaks = peaks, BSgenomeName = Hsapiens)

  testthat::expect_equal(as.character(temp1$aln.all$guideAlignment2OffTarget),
                         "...............^.....")

  testthat::expect_true(as.character(temp1$aln.all$offTarget_sequence) ==
                          as.character(reverseComplement(peaks)))

  testthat::expect_equal(as.numeric(temp1$aln.all$pos.insertion), 16)

  cat("Start testing bulge on gRNA offtarget on minus strand with mismatches...")

  # PAM followed by 2t(19A), 19c(2G), 20t(1A), 6 insertion
  peaks <- DNAStringSet(c("CCTGTGGCTGGGGTGGAGGGGGCT"))
  names(peaks) <- "testMinusBulgeOff"
  temp1 <- getAlnWithBulge(gRNA, gRNA.name = "T2",
                           peaks = peaks, BSgenomeName = Hsapiens)

  testthat::expect_equal(as.character(temp1$aln.all$guideAlignment2OffTarget),
                         "AG.............^...A.")

  testthat::expect_true(as.character(reverseComplement(DNAStringSet(substr(peaks,
           temp1$aln.all$offTarget_Start, temp1$aln.all$offTarget_End)))) ==
                          temp1$aln.all$offTarget_sequence)

  testthat::expect_true(as.character(temp1$aln.all$offTarget_sequence) ==
                          as.character(reverseComplement(peaks)))

  testthat::expect_equal(as.numeric(temp1$aln.all$pos.insertion), 16)

  testthat::expect_equal(as.character(temp1$aln.all$PAM.sequence), "AGG")
  testthat::expect_equal(
    as.numeric(unlist(temp1$aln.all$pos.mismatch)),c(19,2,1))

  #### hg19
  #detach("package:BSgenome.Hsapiens.UCSC.hg38", unload = TRUE)
  require(BSgenome.Hsapiens.UCSC.hg19)
  peaks.f <- system.file("extdata", "T2plus100OffTargets.bed",
                        package = "GUIDEseq")
  gRNA <- substr(as.character(readDNAStringSet(system.file("extdata", "T2.fa",
                      package = "CRISPRseek"))),
                 1, 20)
  # gRNA 1618 from Linda
  #gRNA <- "TTGCTTTTATCACAGGCTCC"
  # peaks with hg38 from Linda
  names(gRNA) <- "T2"
  temp <- getAlnWithBulge(gRNA, gRNA.name = "T2",
                          peaks = peaks.f, BSgenomeName = Hsapiens)

  bed <- read.table(system.file("extdata", "T2plus100OffTargets.bed",
                                package = "GUIDEseq"),
                    sep = "\t",
                    header = TRUE)
  merged.bed <- merge(bed, temp$aln.all, by.x = "names", by.y ="offTarget")
  # some sequences do not contain the PAM.pattern
  #testthat::expect_equal(nrow(merged.bed), nrow(temp))
  testthat::expect_equal(merged.bed$peak_score,merged.bed$totalCount)
  testthat::expect_equal(merged.bed$chr.x,merged.bed$chr.y)

  off.names <- c("chr1+:121373900:121373925 chr1-:121373840:121373865",
                 "chr2+:221535820:221535845 chr2-:221535730:221535755",
                 "chr10+:42529870:42529895 chr10-:42529770:42529795")

  cat("Start testing offtarget on the plus strand ...")
  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names == off.names[1],]$pos.indel))
    , 7)
  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[1],]$guideAlignment2OffTarget)),
    "..G.GG-....AAT.T..A.")
  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[1],]$n.guide.mismatch)), 8)
  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names == off.names[1],]$PAM.sequence)),
    "TTG")
  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[1],]$pos.mismatch)),
    c(3, 5, 6, 12, 13, 14, 16, 19))

  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[1],]$offTarget_End)),
    121373887)

  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[1],]$offTarget_Start)),
    121373866)

  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[1],]$chromosome.y)),
    "chr1")

  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[1],]$mismatch.distance2PAM)),
    "18,16,15,9,8,7,5,2")

  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                     off.names[1],]$n.guide.mismatch)),
    8)

  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[1],]$offTarget_sequence)),
    "GAGCGG-TCCAAATCTCCACTTG")

  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[1],]$offTargetStrand)),
    "+")

  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[3],]$offTarget_Start)),
    42529795)

  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[3],]$offTarget_End)),
    42529816)

  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[3],]$chromosome.y)),
    "chr10")


  cat("Start testing offtarget on the minus strand ...")

  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[2],]$offTarget_sequence)),
    "CACCCCC-CCACCCCACCCCAGG")
  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[2],]$offTargetStrand)),
    "-")

  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                     off.names[2],]$n.guide.mismatch)),
    3)

  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[2],]$mismatch.distance2PAM)),
    "2,5,20")


    testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names == off.names[2],]$pos.indel))
    , 8)

  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                     off.names[2],]$guideAlignment2OffTarget)),
         "C......-.......A..C.")

  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[2],]$n.guide.mismatch)), 3)
  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names == off.names[2],]$PAM.sequence)),
     "AGG")
  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[2],]$pos.mismatch)),
       c(19, 16, 1))



  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[2],]$offTarget_Start)),
    221535826)

  testthat::expect_equal(
    as.numeric(unlist(merged.bed[merged.bed$names ==
                                   off.names[2],]$offTarget_End)),
    221535847)

  testthat::expect_equal(
    as.character(unlist(merged.bed[merged.bed$names ==
                                   off.names[2],]$chromosome.y)),
    "chr2")
})

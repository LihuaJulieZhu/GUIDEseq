detach("package:BSgenome.Hsapiens.UCSC.hg19", unload = TRUE)
test_that("getOfftargetScoreBulge with deletion in gRNA works", {
  library(BSgenome.Hsapiens.UCSC.hg38)
  peaks <- system.file("extdata", "1450-chr14-chr2-bulge-test.bed",
      package = "GUIDEseq")
  mismatch.activity.file <- system.file("extdata", "NatureBiot2016SuppTable19DoenchRoot.xlsx",
         package = "GUIDEseq")
  gRNA <- "TGCTTGGTCGGCACTGATAG"
  gRNA.name <- "Test1450"

  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
  temp <- offTargetAnalysisWithBulge(gRNA = gRNA, gRNA.name = gRNA.name,
        peaks = peaks, BSgenomeName = Hsapiens, mat = mat,
        mismatch.activity.file = mismatch.activity.file)
  expect_equal(c(1.000000, 0.028614), 
    temp$score.bulges$predicted_cleavage_score)
})

detach("package:BSgenome.Hsapiens.UCSC.hg38", unload = TRUE)
library(BSgenome.Hsapiens.UCSC.hg19)

test_that("getOfftargetScoreBulge with insertion in gRNA works", {
    peaks.f <- system.file("extdata", "T2plus100OffTargets.bed",
                       package = "GUIDEseq")
   gRNA <-"GACCCCCTCCACCCCGCCTC"
   gRNA.name  <- "T2"
   temp <- offTargetAnalysisWithBulge(gRNA = gRNA, gRNA.name = gRNA.name,
                      peaks = peaks.f, BSgenomeName = Hsapiens,
                      mismatch.activity.file = mismatch.activity.file,
                      peaks.withHeader = TRUE)
    expect_equal(0.393846 , temp$score.bulges$predicted_cleavage_score[11])
})


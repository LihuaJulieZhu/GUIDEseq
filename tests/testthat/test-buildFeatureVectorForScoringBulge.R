#devtools::load_all("~/Dropbox (UMass Medical School)/Bioconductor/Trunk/GUIDEseq/")
library(BSgenome.Hsapiens.UCSC.hg19)

test_that("buildFeatureVectorsForScoring works", {
  gRNA.file <- system.file("extdata", "gRNA.fa", package = "GUIDEseq")
peaks <- system.file("extdata", "gRNA-PlusMinusPeaksMerged.bed",
                     package = "GUIDEseq")
gRNA <- readDNAStringSet(gRNA.file, use.names = TRUE)
gRNA.name <- names(gRNA)
gRNA <- as.character(gRNA)
#peaks <- "gRNA-PlusMinusPeaksMerged.bed"
gRNA <- substr(gRNA, 1,20)
alns <- getAlnWithBulge(gRNA, gRNA.name = gRNA.name,
                        peaks = peaks, BSgenomeName = Hsapiens)
fv <- buildFeatureVectorForScoringBulge(alns$aln.all)
# important since getOfftargetScore2 expects a colname with name aligment
colnames(fv$featureVectors)[colnames(fv$featureVectors) ==
                                    "guideAlignment2OffTarget"] <- "alignment"
featureVectors <- CRISPRseek:::getOfftargetScore2(fv$featureVectors)
featureVectors <- featureVectors[order(featureVectors$offTarget), ]
featureVectors <- featureVectors[, -(grep("pos", colnames(featureVectors)))]

expect_warning(t2 <- offTargetAnalysisOfPeakRegions(gRNA.file, BSgenomeName = Hsapiens,
                                     peaks = "gRNA-PlusMinusPeaksMerged.bed",
                                     outputDir = "testoff",
                                     scoring.method = "CFDscore")
)
mismatch.activity.file <-
  system.file("extdata","NatureBiot2016SuppTable19DoenchRoot.xlsx",
              package = "GUIDEseq")

res.bulge <- offTargetAnalysisWithBulge(gRNA = gRNA, gRNA.name = gRNA.name,
                                        peaks = "gRNA-PlusMinusPeaksMerged.bed",
                                        BSgenomeName = Hsapiens,
                                        mismatch.activity.file = mismatch.activity.file)

t2 <- t2[order(t2$names), ]

res.bulge <- as.data.frame(res.bulge$score.bulges)
res.bulge <- res.bulge[order(res.bulge$offTarget), ]

res.bulge <- res.bulge[res.bulge$n.insertion == 0 & res.bulge$n.deletion == 0,]

expect_equal(as.numeric(res.bulge$predicted_cleavage_score),
             as.numeric(t2$predicted_cleavage_score))
expect_equal(as.numeric(res.bulge$offTarget_Start),
             as.numeric(t2$offTarget_Start))

expect_equal(as.numeric(res.bulge$offTarget_End),
             as.numeric(t2$offTarget_End))

expect_equal(unlist(res.bulge$offTargetStrand),
             t2$offTargetStrand)

expect_equal(unlist(res.bulge$offTarget_sequence),
             t2$offTarget_sequence)

expect_equal(unlist(res.bulge$chromosome),
             t2$chromosome)
})



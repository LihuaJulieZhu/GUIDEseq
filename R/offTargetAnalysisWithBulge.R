.updateScore <- function(fv.geThan1, fv.lessThan1,
                         col.prefix = "IsInsertion.pos",
                         weights, position, type)
{
    if (dim(fv.geThan1)[1] > 0)
    {
        indel.pos <- fv.geThan1[,grep(col.prefix,
                                      colnames(fv.geThan1))]
        indel.pos <- apply(indel.pos, 1, as.numeric)
        pos <- grep(col.prefix, colnames(fv.geThan1))
        min.pos <- min(pos)

        if (length(pos) > 0)
        {
            score.new <- unlist(lapply(1:nrow(fv.geThan1), function(i)
            {
                indel.index <- pos[fv.geThan1[i, pos] == 1] - min.pos + 1
                if (col.prefix == "IsInsertion.pos")
                    this.type <- unlist(strsplit(as.character(
                        fv.geThan1[i,]$gRNA.insertion), ","))
                else
                    this.type <- unlist(strsplit(as.character(
                        fv.geThan1[i,]$gRNA.deletion), ","))
                score.new1 <- fv.geThan1[i, ]$score
                for (j in 1:length(indel.index))
                {
                    if (indel.index[j] > 1)
                        score.new1 <- score.new1 *
                            weights[position ==
                                indel.index[j] &
                                type == this.type[j]]
                }
                score.new1
            }))
            fv.geThan1$score <- score.new
        }
    }
    if (dim(fv.geThan1)[1] > 0)
    {
        score <- fv.geThan1
        if (dim(fv.lessThan1)[1] > 0)
            score <- rbind(fv.lessThan1, fv.geThan1)
    }
    else
    {
        score <- fv.lessThan1
    }
    score
}

#' @importFrom BiocGenerics subset unlist lapply rbind
#' @importFrom hash hash values
#' @importFrom rio import
#' @importFrom hash hash
#'
#' @example
#'   peaks <- system.file("extdata","1450-chr14-chr2-bulge-test.bed", package = "GUIDEseq")
#'   mismatch.activity.file <-system.file("extdata", "NatureBiot2016SuppTable19DoenchRoot.xlsx", 
#'     package = "GUIDEseq")
#'
#'   gRNA <- "TGCTTGGTCGGCACTGATAG"
#'   gRNA.name <- "Test1450"
#'   library(BSgenome.Hsapiens.UCSC.hg38)
#'
#'   temp <- offTargetAnalysisWithBulge(gRNA = gRNA, gRNA.name = gRNA.name,
#'      peaks = peaks, BSgenomeName = Hsapiens,
#'      mismatch.activity.file = mismatch.activity.file)
#'
offTargetAnalysisWithBulge <-
    function(gRNA, gRNA.name,
             peaks,
             BSgenomeName,
             mat,
             peaks.withHeader = FALSE,
             peaks.format= "bed",
             gapOpening = 1L,
             gapExtension = 3L,
             max.DNA.bulge = 2L,
             max.mismatch = 10L,
             allowed.mismatch.PAM = 2L,
             upstream = 20L,
             downstream =20L,
             PAM.size = 3L,
             gRNA.size = 20L,
             PAM = "NGG",
             PAM.pattern = "NNN$",
             PAM.location = "3prime",
        mismatch.activity.file = system.file("extdata",
            "NatureBiot2016SuppTable19DoenchRoot.xlsx",
            package = "CRISPRseek")
)
{
    alns <- getAlnWithBulge(gRNA, gRNA.name = gRNA.name,
            peaks = peaks, BSgenomeName = BSgenomeName,
            peaks.withHeader = peaks.withHeader,
            peaks.format = peaks.format,
            gapOpening = gapOpening,
            gapExtension = gapExtension,
            max.DNA.bulge = max.DNA.bulge,
            max.mismatch = max.mismatch ,
            allowed.mismatch.PAM = allowed.mismatch.PAM,
            upstream = upstream ,
            downstream = downstream,
            PAM.size = PAM.size,
            gRNA.size = gRNA.size,
            PAM = PAM,
            PAM.pattern = PAM.pattern,
            PAM.location = PAM.location
            )
    temp <- alns$aln.all
    alns <- subset(temp,
                   (as.numeric(temp$n.mismatch) +
                        as.numeric(temp$n.insertion) +
                        as.numeric(temp$n.deletion)) < max.mismatch)

    fv <- buildFeatureVectorForScoringBulge(alns)
    colnames(fv$featureVectors)[colnames(fv$featureVectors) ==
                                    "guideAlignment2OffTarget"] <- "alignment"
    featureVectors <- CRISPRseek:::getOfftargetScore2(fv$featureVectors)

    scores.mismatch <-  featureVectors$score
    deletion.activity <- import(mismatch.activity.file, which = 2)
    insertion.activity <- import(mismatch.activity.file, which = 3)
    required.col.ins <- c("Insertion.Type", "Position", "Percent.Active")
    required.col.del <- c("Deletion.Type", "Position", "Percent.Active")
    if (length(intersect(colnames(insertion.activity), required.col.ins)) !=
        length(required.col.ins))
    {
         stop("Please rename the insertion sheet of the activity file column
         to contain at least these 3 column names: Insertion.Type,
              Position, Percent.Active\n")
    }
    if (length(intersect(colnames(deletion.activity), required.col.del)) !=
        length(required.col.del))
    {
        stop("Please rename the deletion sheet of the activity file column
         to contain at least these 3 column names: Deletion.Type,
              Position, Percent.Active\n")
    }
    colnames(insertion.activity)[colnames(insertion.activity) ==
                                     "Insertion.Type"] <- "Type"

    colnames(deletion.activity)[colnames(deletion.activity) ==
                                     "Deletion.Type"] <- "Type"

    #mismatch.activity[mismatch.activity$Mismatch.Type == "rA:dG" &
    #    mismatch.activity$Position == 10,]$Percent.Active
    ##### by default weights is a column vector
    ##### the mismatch activity  is given as pos 20, 19, 18,....1 distance from PAM,    ##### Position named as 1, 2, 3, .....20 though
    ##### and the featureVectors is in the same order now
    ##### so no need to reverse any more. weights = rev(weights)
    fv.geThan1 <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.insertion)) >= 1)
    fv.lessThan1 <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.insertion)) < 1)

    featureVectors <- .updateScore(fv.geThan1, fv.lessThan1,
                          col.prefix = "IsInsertion.pos",
                          weights = insertion.activity$Percent.Active,
                          type = insertion.activity$Type,
                          position = insertion.activity$Position)

    # do the same for deletion
    fv.geThan1 <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.deletion)) >= 1)
    fv.lessThan1 <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.deletion)) < 1)

    score <- .updateScore(fv.geThan1, fv.lessThan1,
                    col.prefix = "IsDeletion.pos",
                    weights = deletion.activity$Percent.Active,
                    type = deletion.activity$Type,
                    position = deletion.activity$Position)

    score$alignment <- as.character(score$alignment)
    score$score <- round(score$score, 6)
    colnames(score)[colnames(score) == "score"] <- "predicted_cleavage_score"
    score  <- score[order(c(score$name,score$predicted_cleavage_score),
                          decreasing = TRUE), ]
    score <- score[, -grep(".pos", colnames(score))]
    colnames(score)[colnames(score) == "alignment"] <-
                                    "guideAlignment2OffTarget"
    list(score.bulges = unique(score[!is.na(score$predicted_cleavage_score), ]),
         offtargets.dropped = fv$featureVectors.dropped,
         scores.mismatch =  scores.mismatch )
}

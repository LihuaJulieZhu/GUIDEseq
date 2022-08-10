#' @title Compare Samples using Fisher's exact test
#' @param df a data frame containing the peak score and sequence depth for each
#' sample
#' @param col.count1 the score (e.g., peak_score) column used as the numerator
#' for calculating odds ratio. For example,if the tenth column contains the score
#' for sample 1, then set col.count1 = 10
#' @param col.count2 the score (e.g., peak_score) column used as the denominator
#' for calculating odds ratio. For example,if the nineteenth column contains the score
#' for sample 1, then set col.count2 = 19
#' @param total1 the sequence depth for sample 1
#' @param total2 the sequence depth for sample 2
#' @param multiAdjMethod A vector of character strings containing the names of
#' the multiple testing procedures for which adjusted p-values are to be
#' computed. This vector should include any of the following:
#' "none", "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY",
#' "ABH", and "TSBH". Please type ?multtest::mt.rawp2adjp for details.
#' Default to "BH"
#' @param comparison.score the score to be used for statistical analysis.
#' Two options are available: "peak_score" and "umi.count"
#' umi.count is the number of unique UMIs in the associated peak region
#' without considering the sequence coordinates while peak_score takes
#' into consideration of the sequence coordinates
#'
#' @author Lihua Julie Zhu
#' @importFrom multtest mt.rawp2adjp
#' @importFrom stringr str_c
#' @importFrom methods is
#' @importFrom stats fisher.test
#'

compareSamples <- function(df,
                           col.count1,
                           col.count2,
                           total1,
                           total2,
                           multiAdjMethod = "BH",
                           comparison.score = c("peak_score", "umi.count"))
{
   stopifnot(!missing(df),
             class(df) == "data.frame",
             nrow(df) >=1,
             !missing(col.count1),
             !missing(col.count2),
             !missing(total1),
             !missing(total2))
  comparison.score <- match.arg(comparison.score)

   df[is.na(df[, col.count1]), col.count1] <- 0
   df[is.na(df[, col.count2]), col.count2] <- 0

   group1 <- colnames(df)[col.count1]
   group2 <- colnames(df)[col.count2]
   sample1 <- gsub(paste0(".", comparison.score), "", group1, fixed = TRUE)
   sample2 <- gsub(paste0(".", comparison.score), "", group2, fixed = TRUE)

   unique.pair <- unique(cbind(as.numeric(df[, col.count1]),
                           as.numeric(df[, col.count2])))
   colnames(unique.pair) = c(group1, group2)
   tryCatch(pvalue <- do.call(rbind, lapply(1:dim(unique.pair)[1], function(i) {
      temp = matrix(c(as.numeric(unique.pair[i, 1]),
                    as.numeric(unique.pair[i, 2]),
                  total1 - as.numeric(unique.pair[i, 1]),
                  total2 - as.numeric(unique.pair[i, 2])), nrow = 2,
                  dimnames = list(c(group1, group2), c("yes", "no")))
      r1 = fisher.test(temp)
      c(unique.pair[i, ], r1$p.value, r1$est)
      })), error = function(e) {print("Some peak score is larger than total1 or total2 is too small. Please double check!")})

   colnames(pvalue) = c(group1,
                     group2,
                     "p.value",
                      str_c(c(sample1, "vs", sample2,"odds.ratio"),
                                collapse = "."))
   df <- merge(df, pvalue)
   if (multiAdjMethod != "none") {
        procs = c(multiAdjMethod)
        ptemp = df[, dim(df)[2] - 1]
        res <- mt.rawp2adjp(ptemp, procs)
        adjp = unique(res$adjp)
        colnames(adjp)[1] =  "p.value"
        colnames(adjp)[2] = str_c(c(sample1, "vs", sample2, multiAdjMethod,
                              "adjusted.p.value"),
                            collapse  = ".")
        df = merge(df, adjp, all.x = TRUE)
        colnames(df)[1] = str_c(c(sample1, "vs", sample2, "p.value"),
                                collapse =".")
        df[is.na(df[, dim(df)[2]]),dim(df)[2]] = 1
        df <- cbind(df[, -c(1:3)], df[,1:3])
  }
  else {
        colnames(df)[dim(df)[2] - 1] <- str_c(c(group1, "vs",
                                                        group2, "p.value"),
                                             collapse =".")
        df <- cbind(df[, -c(1:2)], df[,1:2])
  }
  df
}

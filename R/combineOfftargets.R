#' Combine Offtargets
#'
#' Merge offtargets from different samples
#'
#' Please note that by default, merged file will only contain peaks with
#' offtargets found in the genome in GUIDEseqAnalysis function.
#'
#' @param offtarget.folder offtarget summary output folders created in
#' GUIDEseqAnalysis function
#' @param sample.name Sample names to be used as part of the column names in
#' the final output file
#' @param remove.common.offtargets Default to FALSE If set to TRUE, off-targets
#' common to all samples will be removed.
#' @param control.sample.name The name of the control sample for filtering
#' off-targets present in the control sample
#' @param offtarget.filename Default to offTargetsInPeakRegions.xls, generated
#' in GUIDEseqAnalysis function
#' @param common.col common column names used for merge files. Default to
#' c("total.mismatch.bulge","chromosome", "offTarget_Start","offTarget_End",
#' "offTargetStrand","offTarget_sequence","PAM.sequence","guideAlignment2OffTarget",
#' "mismatch.distance2PAM","n.guide.mismatch","n.PAM.mismatch",
#' "n.DNA.bulge","n.RNA.bulge","pos.DNA.bulge","DNA.bulge","pos.RNA.bulge",
#' "RNA.bulge","gRNA.name","gRNAPlusPAM","predicted_cleavage_score",
#' "inExon","symbol","entrez_id")
#' @param exclude.col columns to be excluded before merging. Please check
#' offTargetsInPeakRegions.xls to choose the desired columns to exclude
#' @param outputFileName The merged offtarget file
#' @param comparison.sample1 A vector of sample names to be used for comparison.
#' For example, comparison.sample1 = c("A", "B"),
#' comparison.sample2 = rep("Control", 2) indicates that you are
#' interested in comparing sample A vs Control and B vs Control
#' Please make sure the sample names specified in comparison.sample1 and
#' comparison.sample2 are in the sample name list specified in sample.name
#' @param comparison.sample2 A vector of sample names to be used for comparison.
#'For example, comparison.sample1 = c("A", "B"),
#' comparison.sample2 = rep("Control", 2) indicates that you are
#' interested in comparing sample A vs Control and B vs Control
#' @param multiAdjMethod A vector of character strings containing the names of
#' the multiple testing procedures for which adjusted p-values are to be
#' computed. This vector should include any of the following: "none",
#' "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH",
#' and "TSBH". Please type ?multtest::mt.rawp2adjp for details. Default to "BH"
#' @param comparison.score the score to be used for statistical analysis.
#' Two options are available: "peak_score" and "umi.count"
#' umi.count is the number of unique UMIs in the associated peak region
#' without considering the sequence coordinates while peak_score takes
#' into consideration of the sequence coordinates
#' @param overwrite Indicates whether to overwrite the existing file
#' specified by outputFileName, default to FALSE.
#' @return a data frame containing all off-targets from all samples merged by
#' the columns specified in common.col. Sample specific columns have sample.name
#' concatenated to the original column name, e.g., peak_score becomes
#' sample1.peak_score.
#' @author Lihua Julie Zhu
#' @keywords misc
#' @examples
#'
#'     offtarget.folder <- system.file("extdata",
#'         c("sample1-17", "sample2-18", "sample3-19"),
#'         package = "GUIDEseq")
#'     mergedOfftargets <-
#'        combineOfftargets(offtarget.folder = offtarget.folder,
#'        sample.name = c("Cas9Only", "WT-SpCas9", "SpCas9-MT3-ZFP"),
#'	  comparison.sample1 = c("Cas9Only", "SpCas9-MT3-ZFP"),
#'	  comparison.sample2 = rep("WT-SpCas9", 2),
#'        outputFileName = "TS2offtargets3Constructs.xlsx")
#'
#' @importFrom utils read.table write.table
#' @importFrom limma vennCounts vennDiagram
#' @importFrom openxlsx createWorkbook saveWorkbook addWorksheet createStyle writeData
#'
#' @export combineOfftargets
#'
combineOfftargets <- function(offtarget.folder,
    sample.name, remove.common.offtargets = FALSE,
    control.sample.name,
    offtarget.filename = "offTargetsInPeakRegions.xls",
    common.col = c("total.mismatch.bulge",
                   "chromosome",
                   "offTarget_Start",
                   "offTarget_End",
                   "offTargetStrand",
                   "offTarget_sequence",
                   "PAM.sequence",
                   "guideAlignment2OffTarget",
                   "mismatch.distance2PAM",
                   "n.guide.mismatch",
                   "n.PAM.mismatch",
                   "n.DNA.bulge",
                   "n.RNA.bulge",
                   "pos.DNA.bulge",
                   "DNA.bulge",
                   "pos.RNA.bulge",
                   "RNA.bulge",
                   "gRNA.name",
                   "gRNAPlusPAM",
                   "predicted_cleavage_score",
                   "inExon",
                   "symbol",
                   "entrez_id"),
    exclude.col = "",
    outputFileName,
    comparison.sample1,
    comparison.sample2,
    multiAdjMethod = "BH",
    comparison.score = c("peak_score", "umi.count"),
    overwrite = FALSE)
{
    stopifnot(!missing(offtarget.folder), length(offtarget.folder) > 1,
              !missing(sample.name), !missing(outputFileName))
     if(!missing(comparison.sample1) &&
          (length(setdiff(comparison.sample1, sample.name)) > 0 ||
           missing(comparison.sample2) ||
          length(comparison.sample1) != length(comparison.sample2) ||
          length(setdiff(comparison.sample2, sample.name)) > 0 ))
        stop("Please make sure that comparison.sample1 has the same number of
             samples as that of comparison.sample2 with the names
             listed in the sample.name!")
    comparison.score <- match.arg(comparison.score)
    for (i in 1:length(offtarget.folder)) {
        if(!file.exists(file.path(offtarget.folder[i],
                    offtarget.filename, fsep = .Platform$file.sep)))
            stop(offtarget.filename, " is not found in the directory ",
                 offtarget.folder[i])
    }
    if(missing(outputFileName) ||
       (file.exists(outputFileName) && !overwrite))
        stop("outputFileName is required and
              cannot already exists if overwrite is set to FALSE!")

    all <- read.table(file.path(offtarget.folder[1],
        offtarget.filename, fsep = .Platform$file.sep),
        sep="\t", header = TRUE, stringsAsFactors = FALSE)
    all <- subset(all, !is.na(all$offTarget))
    if (!missing(exclude.col) && exclude.col != "")
        all <- all[, -which(colnames(all) %in% exclude.col)]
    # if(length(setdiff(common.col, colnames(all))) > 0)
    # {
    #     message(paste(paste(setdiff(common.col, colnames(all)),collapse = ", "),
    #         "are not valid column names! Please make sure to
    #         specify the correct column names for common.col! If
    #         you have specified the column names in exclued.col, please
    #         do not specify in the common.col!"))
    # }

    common.col <- intersect(setdiff(common.col, c("offTarget", "peak_score")),
                      colnames(all))
    colnames(all)[!colnames(all) %in% common.col] <- paste(sample.name[1],
        colnames(all)[!colnames(all) %in% common.col], sep=".")

    if (length(grep("pos.RNA.bulge", colnames(all))) >0 )
    {
        all$pos.RNA.bulge[is.na(all$pos.RNA.bulge)] <- ""
        all$pos.DNA.bulge[is.na(all$pos.DNA.bulge)] <- ""
    }
    for (i in 2:length(offtarget.folder))
    {
        off <- read.table(file.path(offtarget.folder[i], offtarget.filename,
            fsep = .Platform$file.sep), sep="\t", header = TRUE,
            stringsAsFactors = FALSE)
        off <- subset(off, !is.na(off$offTarget_sequence))
        if (!missing(exclude.col) && exclude.col != "")
            off <- off[,-which(colnames(off) %in% exclude.col)]
        colnames(off)[!colnames(off) %in% common.col] <- paste(sample.name[i],
            colnames(off)[!colnames(off) %in% common.col], sep=".")
        if (length(grep("pos.RNA.bulge", colnames(off))) > 0 )
        {
           off$pos.RNA.bulge[is.na(off$pos.RNA.bulge)] <- ""
           off$pos.DNA.bulge[is.na(off$pos.DNA.bulge)] <- ""
        }
        all <- merge(all, off, by = common.col, all = TRUE)
    }

    offtarget.columns <- paste(sample.name, "peak_score",
                               sep = "." )
    temp <- do.call(cbind,
        lapply(offtarget.columns, function(i) {
            temp1 <- all[, i]
            as.numeric(!is.na(temp1))
        }))
    colnames(temp) <- sample.name
    venn_cnt <- vennCounts(temp)

    vennDiagram(venn_cnt)
    if (remove.common.offtargets)
    {
        all <- subset(all, rowSums(temp) < dim(temp)[2])
    }
    if (!missing(control.sample.name))
    {
        if (control.sample.name %in%  sample.name)
            all <- subset(all, !temp[, control.sample.name])
        else
            message("Please note that control.sample.name is not on the sample.name list, filtering skipped!")
    }
    #write.table(subset(all, !is.na(all[,1])), file = outputFileName, sep="\t",
      #          row.names=FALSE)

    if (length(grep("sequence.depth", colnames(all))) > 0)
    {
        for (i in 1:length(comparison.sample1))
        {
            # important to have the following 4 lines inside the loop because
            # the order of the columns may change after calling compareSamples
            col.count1 <- which(colnames(all) ==
                                paste(comparison.sample1[i], comparison.score,
                                      sep = "."))
            col.count2 <- which(colnames(all) ==
                                paste(comparison.sample2[i], comparison.score,
                                      sep = "."))
            col.total1 <-  which(colnames(all) %in%
                                 paste(comparison.sample1[i],"sequence.depth",
                                       sep = "."))
            col.total2 <- which(colnames(all) %in%
                                paste(comparison.sample2[i],"sequence.depth",
                                      sep = "."))
            total1 = max(all[, col.total1], na.rm = TRUE)
            total2 = max(all[, col.total2], na.rm = TRUE)

            all[is.na(all[, col.total1]),col.total1] <- total1
            all[is.na(all[, col.total2]),col.total2] <- total2

            all <- compareSamples(all, col.count1 = col.count1,
                       col.count2 = col.count2,
                       total1 = total1,
                       total2 = total2,
                       multiAdjMethod = multiAdjMethod,
                       comparison.score = comparison.score)
        }
        # write.table(all, file = outputFileName,
        #     sep="\t", row.names=FALSE)
    }
    workbook <- createWorkbook()
    addWorksheet(workbook, "mergedOfftargets")
    #bold.style <- createStyle(textDecoration = c("underline","Bold"))
    bold.style <- createStyle(textDecoration = "Bold")
    writeData(workbook, "mergedOfftargets", all, startRow = 1, startCol = 1,
              headerStyle = bold.style)
    saveWorkbook(workbook, file = outputFileName,
                 overwrite = overwrite)
    all
}

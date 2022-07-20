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
#' c("offTarget","predicted_cleavage_score", "gRNA.name", "gRNAPlusPAM",
#' "offTarget_sequence", "guideAlignment2OffTarget", "offTargetStrand",
#' "mismatch.distance2PAM", "n.PAM.mismatch", "n.guide.mismatch",
#' "PAM.sequence", "offTarget_Start", "offTarget_End", "chromosome")
#' @param exclude.col columns to be excluded before merging.  Please check
#' offTargetsInPeakRegions.xls to choose the desired columns to exclude
#' @param outputFileName The merged offtarget file
#' @return a tab-delimited file similar to offTargetsInPeakRegions.tsv,
#' containing all peaks from all samples merged by potential gRNA binding
#' sites, mismatch number and positions, alignment to the input gRNA and
#' predicted cleavage score. Sample specific columns have sample.name
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
#'        sample.name = c("cas9Only", "WT SpCas9", "SpCas9-MT3-ZFP"),
#'        outputFileName = "TS2offtargets3Constructs.xls")
#'
#' @importFrom utils read.table write.table
#' @importFrom limma vennCounts vennDiagram
#'
#' @export combineOfftargets
#'
combineOfftargets <- function(offtarget.folder,
    sample.name, remove.common.offtargets = FALSE,
    control.sample.name,
    offtarget.filename = "offTargetsInPeakRegions.xls",
    common.col = c("offTarget","predicted_cleavage_score",
        "gRNA.name", "gRNAPlusPAM", "offTarget_sequence",
        "guideAlignment2OffTarget", "offTargetStrand",
        "mismatch.distance2PAM", "n.PAM.mismatch",
        "n.guide.mismatch", "PAM.sequence", "offTarget_Start",
        "offTarget_End", "chromosome"),
    exclude.col,
    outputFileName)
{
    all <- read.table(file.path(offtarget.folder[1],
        offtarget.filename, fsep = .Platform$file.sep),
        sep="\t", header = TRUE)
    all <- subset(all, !is.na(all$offTarget))
    if (!missing(exclude.col) && exclude.col != "")
        all <- all[, -which(colnames(all) %in% exclude.col)]
    if(length(setdiff(common.col, colnames(all))) > 0)
    {
        stop(paste(setdiff(common.col, colnames(all)),
            "are not valid column names! Please make sure to
            specify the correct column names for common.col! If
            you have specified the column names in exclued.col, please
            do not specify in the common.col!"))
    }

    common.col <- setdiff(common.col, "peak_score")
    colnames(all)[!colnames(all) %in% common.col] <- paste(sample.name[1],
        colnames(all)[!colnames(all) %in% common.col], sep=".")

    for (i in 2:length(offtarget.folder))
    {
        off <- read.table(file.path(offtarget.folder[i], offtarget.filename,
            fsep = .Platform$file.sep), sep="\t", header = TRUE)
        off <- subset(off, !is.na(off$offTarget))
        if (!missing(exclude.col) && exclude.col != "")
            off <- off[,-which(colnames(off) %in% exclude.col)]
        colnames(off)[!colnames(off) %in% common.col] <- paste(sample.name[i],
            colnames(off)[!colnames(off) %in% common.col], sep=".")

        all <- merge(all, off, by = common.col, all = TRUE)
    }

    offtarget.columns <- paste(sample.name, "peak_score", sep = "." )
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
            warning("Please note that control.sample.name is not on the sample.name list, filtering skipped!")
    }
    write.table(subset(all, !is.na(all[,1])), file = outputFileName, sep="\t", row.names=FALSE)

    all
}

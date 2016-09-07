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
        all <- all[, -which(colnames(all) == exclude.col)]
    if(length(setdiff(common.col, colnames(all))) > 0)
    {
        stop(paste(setdiff(common.col, colnames(all)), 
            "are not valid column names! Please make sure to 
            specify the correct column names for common.col!"))
    }

    colnames(all)[!colnames(all) %in% common.col] <- paste(sample.name[1], 
        colnames(all)[!colnames(all) %in% common.col], sep=".")

    for (i in 2:length(offtarget.folder))
    {
        off <- read.table(file.path(offtarget.folder[i], offtarget.filename,
            fsep = .Platform$file.sep), sep="\t", header = TRUE)
        off <- subset(off, !is.na(off$offTarget))
        if (!missing(exclude.col) && exclude.col != "")
            off <- off[,-which(colnames(off) == exclude.col)]
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

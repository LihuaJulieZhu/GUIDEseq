combineOfftargets <- function(offtarget.folder, 
    sample.name, 
    offtarget.filename = "offTargetsInPeakRegions.xls",
    common.col = c("targetSeqName", "chromosome", "offTargetStrand",
        "offTarget_Start", "offTarget_End","gRNAPlusPAM", 
        "offTarget_sequence", "n.mismatch",  
        "guideAlignment2OffTarget","predicted_cleavage_score"),
    exclude.col = "name",
    outputFileName)
{
    all <- read.table(file.path(offtarget.folder[1], 
        offtarget.filename, fsep = .Platform$file.sep),
        sep="\t", header = TRUE)
    all <- subset(all, !is.na(all$offTarget))
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
        off <- off[,-which(colnames(off) == exclude.col)]
        colnames(off)[!colnames(off) %in% common.col] <- paste(sample.name[i], 
            colnames(off)[!colnames(off) %in% common.col], sep=".")

        all <- merge(all, off, by = common.col, all = TRUE)
    }

    write.table(all, file = outputFileName, sep="\t", row.names=FALSE)
    all
}

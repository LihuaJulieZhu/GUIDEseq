offTargetAnalysisOfPeakRegions <-
    function(gRNA, peaks,
    format=c("fasta", "bed"),
    peaks.withHeader = FALSE,
    BSgenomeName,
    overlap.gRNA.positions = c(17,18),
    upstream = 50,
    downstream =50,
    PAM.size = 3,
    gRNA.size = 20,
    PAM = "NGG",
    PAM.pattern = "(NAG|NGG|NGA)$",
    max.mismatch = 6,
    outputDir,
    allowed.mismatch.PAM = 2,
    overwrite = TRUE,
    weights = c(0, 0, 0.014, 0, 0, 0.395,
    0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,
    0.804, 0.685, 0.583)
   )
{
    thePeaks <- read.table(peaks, sep="\t", header = peaks.withHeader,
        stringsAsFactors = FALSE)
    if (format[2] == "bed")
    {
       if (dim(thePeaks)[2] >= 4 )
           colnames(thePeaks)[1:4] <- c("chromosome", "peak_start", "peak_end",
                "names")
        if(dim(thePeaks)[2] >= 5)
            colnames(thePeaks)[5] = "peak_score"
        if(dim(thePeaks)[2] >= 6)
            colnames(thePeaks)[6] = "peak_strand"
    }
    else
    {
        stop("only bed file with at least 4 columns are supported")
    }
    thePeaks[, 4] <- gsub(" ", "", thePeaks[,4])
    TS2 <- tryCatch(
        (compare2Sequences(inputFile1Path = gRNA, inputFile2Path = peaks,
        findgRNAsWithREcutOnly = FALSE, format = format,
        header = peaks.withHeader,
        BSgenomeName = BSgenomeName,
        upstream = upstream,
        downstream = downstream,
        minREpatternSize = 4,
        findgRNAs = c(FALSE, FALSE),
        removegRNADetails = c(TRUE, FALSE),
        exportAllgRNAs =  "no",
        annotatePaired = FALSE,
        searchDirection = "1to2",
        overlap.gRNA.positions = overlap.gRNA.positions,
        findPairedgRNAOnly = FALSE,
        min.gap = 0, max.gap = 20, gRNA.name.prefix = "gRNA",
        PAM.size = PAM.size,
        gRNA.size = gRNA.size, PAM = PAM, PAM.pattern = PAM.pattern,
        max.mismatch = max.mismatch,
        outputDir = outputDir, foldgRNAs = FALSE,
        allowed.mismatch.PAM = allowed.mismatch.PAM, overwrite = overwrite,
        weights = weights)),
        error = function(e) { print(e); "error" })

    if ( TS2 != "error")
    {
        n.cores <- detectCores() - 1
        cl <- makeCluster(n.cores)
   
        names <- as.character(unlist(parLapply(cl, as.character(TS2$offTarget),
           function(temp) {
           temp1 <- strsplit(temp,":")[[1]]
           end.ind <- length(temp1) - 1
           paste(temp1[1:end.ind], collapse=":")
        })))
        TS2 <- cbind(names = names, TS2)
        excluding.columns = which(colnames(TS2) %in% 
            c("scoreForSeq1", "targetInSeq1", "gRNAefficacy", "scoreDiff"))
        TS2 <- TS2[, -excluding.columns]
        colnames(TS2)[colnames(TS2) == "scoreForSeq2"] = "predicted_cleavage_score"
        colnames(TS2)[colnames(TS2) == "targetInSeq2"] = "offTarget_sequence"
        offTargetOffset <- do.call(rbind, parLapply(cl, as.character(TS2$offTarget),
           function(temp) 
           { 
            temp1 <- strsplit(temp,":")[[1]]
            start.ind <- length(temp1)
            as.numeric(strsplit(temp1[start.ind], "-")[[1]][1:2])
        }))
        stopCluster(cl)
        TS2 <- cbind(TS2, offTarget_Start = offTargetOffset[,1], 
           offTarget_End = offTargetOffset[,2])
        offtargets <- merge(TS2, thePeaks, by = "names", all = TRUE)
        offtargets$offTarget_Start <- 
            ifelse(as.character(offtargets$peak_strand) == "+",
            offtargets$offTarget_Start - upstream + offtargets$peak_start - 1,
            upstream - offtargets$offTarget_End + offtargets$peak_end + 1)
        offtargets$offTarget_End <- 
            offtargets$offTarget_Start + PAM.size  + gRNA.size  - 1
   
        offtargets.minus.minus <- 
            subset(offtargets, as.character(offtargets$peak_strand) == "-" & 
            as.character(offtargets$offTargetStrand) == "-")
        offtargets.minus.plus <- subset(offtargets,
            as.character(offtargets$peak_strand) == "-" & 
            as.character(offtargets$offTargetStrand) == "+")
        if (dim(offtargets.minus.minus)[1] > 0)
            offtargets.minus.minus$offTargetStrand <- "+"
        if (dim(offtargets.minus.plus)[1] > 0)
            offtargets.minus.plus$offTargetStrand<- "-"
   
        offtargets <- rbind(subset(offtargets,
            as.character(offtargets$peak_strand) == "+" | 
            as.character(offtargets$peak_strand) == "*" | 
            is.na(offtargets$offTargetStrand)),
            offtargets.minus.minus, offtargets.minus.plus)
   
        write.table(offtargets, 
            file = file.path(outputDir,"offTargetsInPeakRegions.xls"),
            sep="\t", row.names = FALSE)
        unlink(file.path(outputDir, "scoresFor2InputSequences.xls"))
        offtargets
     }
    else
    {
        write.table(thePeaks,
            file = file.path(outputDir,"offTargetsInPeakRegions.xls"),
            sep="\t", row.names = FALSE)
        thePeaks 
    }
}

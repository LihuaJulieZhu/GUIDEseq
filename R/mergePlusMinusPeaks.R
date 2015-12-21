mergePlusMinusPeaks <-
    function(peaks.gr, peak.height.mcol ="count",
    bg.height.mcol = "bg", distance.threshold = 40, step = 20,
    plus.strand.start.gt.minus.strand.end = TRUE, output.bedfile)
{
    if (missing(peaks.gr)) {
        stop("Missing required peak input peaks.gr! \n")
    }
    if (class(peaks.gr) != "GRanges" ) {
        stop("No valid peaks.gr passed in. It needs to be GRanges object \n")
    }
    if (missing(output.bedfile))
    {
        stop("Missing required bed output file \n")
    } 
    if (length(intersect(names(mcols(peaks.gr)), peak.height.mcol))  == 0)
        stop(paste(peak.height.mcol, "is not a valid metadata column.\n Please",
            "specify a valid metadata column with peak height in peaks.gr \n"))
    names(peaks.gr) <- paste(paste(seqnames(peaks.gr), strand(peaks.gr), sep=""),
        start(peaks.gr),end(peaks.gr), sep=":")
    pos.gr <- subset(peaks.gr, strand(peaks.gr) %in% c( "+", "*"))
    neg.gr <- subset(peaks.gr, strand(peaks.gr) == "-")
    ### peaks from both strand or present in both library
    #if (length(pos.gr) == 0 || length(neg.gr) == 0)
   #     stop("Need peaks from both + and - strand to merge")
    if (length(pos.gr) == 0 || length(neg.gr) == 0 )
    {
	mergedPeaks <- .annotate(from.gr = peaks.gr, to.gr = peaks.gr,
            peak.height.mcol = peak.height.mcol, 
            bg.height.mcol = bg.height.mcol,
            distance.threshold = distance.threshold, step = step,
            plus.strand.start.gt.minus.strand.end = FALSE,
            to.strand = "-")
    }
    else
    {
        mergedPeaks <- .annotate(from.gr = pos.gr, to.gr = neg.gr,
           peak.height.mcol = peak.height.mcol, 
           bg.height.mcol = bg.height.mcol,
           distance.threshold = distance.threshold, step = step,
           plus.strand.start.gt.minus.strand.end = 
           plus.strand.start.gt.minus.strand.end,
           to.strand = "-")
    }
    mergedPeaks.gr <- mergedPeaks$mergedPeaks
    bed <- mergedPeaks$bed
    #mergedPeaksFromPlus2Minus <- mergedPeaks$all.mergedPeaks
    ann.peaks <- mergedPeaks$detailed.mergedPeaks
    if (length(pos.gr) == 0 || length(neg.gr) == 0)
    {
    	peaks.1strandOnly <- peaks.gr[!names(peaks.gr) %in% 
             as.character(ann.peaks$feature) & 
             ! names(peaks.gr) %in% ann.peaks$peak]
    }
    else
    {
        neg.gr2 <- neg.gr[!names(neg.gr) %in% 
            as.character(ann.peaks$feature)]
        pos.gr2 <- pos.gr[!names(pos.gr) %in% ann.peaks$peak]
        peaks.1strandOnly <- c(pos.gr2, neg.gr2)
    }
    write.table(bed, file = output.bedfile, sep = "\t",
        col.names = FALSE, row.names=FALSE, quote=FALSE)
    list(mergedPeaks.gr = mergedPeaks.gr, mergedPeaks.bed = bed,
        peaks.1strandOnly = peaks.1strandOnly, mergedPeaks.details = ann.peaks)
}

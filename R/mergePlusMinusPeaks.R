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
    pos.gr <- subset(peaks.gr, strand(peaks.gr) %in% c( "+", "*"))
    neg.gr <- subset(peaks.gr, strand(peaks.gr) == "-")
    if (length(pos.gr) == 0 || length(neg.gr) == 0)
        stop("Need peaks from both + and - strand to merge")
    names(pos.gr) <- paste(paste(seqnames(pos.gr), strand(pos.gr), sep=""),
        start(pos.gr),end(pos.gr), sep=":")
    names(neg.gr) <- paste(paste(seqnames(neg.gr), strand(neg.gr), sep=""),
        start(neg.gr),end(neg.gr), sep=":")
    mergedPeaks <- .annotate(from.gr = pos.gr, to.gr = neg.gr,
        peak.height.mcol = peak.height.mcol, bg.height.mcol = bg.height.mcol,
        distance.threshold = distance.threshold, step = step,
        plus.strand.start.gt.minus.strand.end = 
        plus.strand.start.gt.minus.strand.end,
        to.strand = "-")
    mergedPeaks.gr <- mergedPeaks$mergedPeaks
    bed <- mergedPeaks$bed
    mergedPeaksFromMinus2Plus = ""
    mergedPeaksFromPlus2Minus <- mergedPeaks$all.mergedPeaks
    if (length(mergedPeaks.gr) > 0)
    {
        ann.peaks <- mergedPeaks$detailed.mergedPeaks
        neg.gr2 <- neg.gr[! names(neg.gr) %in% ann.peaks$feature]
        pos.gr2 <- pos.gr[ names(pos.gr) %in% ann.peaks$peak]
        if (length(neg.gr2) > 0 ) {
            mergedPeaks2 <- .annotate(from.gr = neg.gr2, to.gr = pos.gr2,
                peak.height.mcol = peak.height.mcol,
                bg.height.mcol = bg.height.mcol,
                distance.threshold = distance.threshold, step = step,
                plus.strand.start.gt.minus.strand.end =
                plus.strand.start.gt.minus.strand.end,
                to.strand = "+")
            if (length(mergedPeaks2) > 0)
            {
                bed <- rbind(mergedPeaks$bed, mergedPeaks2$bed)
                mergedPeaks.gr <- unique(c(mergedPeaks$mergedPeaks,
                    mergedPeaks2$mergedPeaks))
                mergedPeaksFromMinus2Plus = mergedPeaks2$all.mergedPeaks
            }
        }
        write.table(bed, file = output.bedfile, sep = "\t",
            col.names = FALSE, row.names=FALSE, quote=FALSE)
        list(mergedPeaks.gr = mergedPeaks.gr, mergedPeaks.bed = bed,
            mergedPeaksFromPlus2Minus = mergedPeaksFromPlus2Minus,
            mergedPeaksFromMinus2Plus = mergedPeaksFromMinus2Plus)
   }
}

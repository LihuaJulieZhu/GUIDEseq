.annotate <-
    function(from.gr, to.gr, peak.height.mcol ="count",
    bg.height.mcol = "bg", distance.threshold = 40, step = 20,
    plus.strand.start.gt.minus.strand.end = TRUE, to.strand = "-")
{
    gr <- annotatePeakInBatch(from.gr, featureType = "TSS", 
        AnnotationData = to.gr, output="nearestStart",
        PeakLocForDistance = "middle", FeatureLocForDistance = "middle")
    gr <- subset(gr, !is.na(gr$distancetoFeature))
    if (plus.strand.start.gt.minus.strand.end) 
        ann.peaks <- as.data.frame(gr[abs(gr$distancetoFeature) <= 
            distance.threshold & gr$distancetoFeature <= step,])
    else
        ann.peaks <- as.data.frame(gr[abs(gr$distancetoFeature) <=
            distance.threshold, ])
    if (dim(as.data.frame(ann.peaks))[1] > 0)
    {
        to.peaks <- as.data.frame(to.gr)
        to.peaks <- cbind(feature = names(to.gr), to.peaks)
        metadata.col <- which(colnames(to.peaks) %in% names(mcols(to.gr)))
        colnames(to.peaks)[metadata.col] <-
            paste(to.strand, colnames(to.peaks)[metadata.col], sep=":")
        to.peaks <- to.peaks[, c(1, metadata.col)]
        temp1 <- merge(to.peaks, ann.peaks)
        temp1 <- cbind(temp1, 
            totalCount = rowSums(temp1[,grep(peak.height.mcol,
            colnames(temp1))]))
        temp1$names <- paste(temp1$peak, temp1$feature, sep=":")
        temp1$minStart <- rowMins(as.matrix(
            temp1[, c("start_position", "start")]))
        temp1$maxEnd <- rowMaxs(as.matrix(temp1[, c("end_position", "end")]))
        bed.temp <- temp1[, c("seqnames", "minStart",
            "maxEnd", "names", "totalCount")]
        bed.temp <- cbind(bed.temp, strand = "+")
        if (length(intersect(names(mcols(to.gr)), bg.height.mcol))  > 0)
        {
            if(is.null(dim(temp1)))
            {
                temp1 <- c(temp1, sum(temp1[,grep(bg.height.mcol,
                    colnames(temp1))]))
            }
            else
            {
                temp1 <- cbind(temp1,
                    totalBg = rowSums(
                    temp1[,grep(bg.height.mcol, colnames(temp1))]))
            }
            mergedPeaks.gr <- GRanges(IRanges(start = 
                as.numeric(as.character(bed.temp[,2])),
                end = as.numeric(as.character(bed.temp[,3])),
                names = bed.temp[,4]),
                seqnames = bed.temp[,1], strand = Rle("+", dim(bed.temp)[1]),
                count = as.numeric(as.character(bed.temp[,5])),
                bg = as.numeric(as.character(temp1$totalBg)))
        }
        else
        {
            mergedPeaks.gr <- GRanges(IRanges(
                start = as.numeric(as.character(bed.temp[,2])),
                end = as.numeric(as.character(bed.temp[,3])),
                names = bed.temp[,4]),
                seqnames = bed.temp[,1], strand = Rle("+", dim(bed.temp)[1]),
                count = as.numeric(as.character(bed.temp[,5])))
        }
        return(list(mergedPeaks = mergedPeaks.gr, bed = bed.temp,
            detailed.mergedPeaks = temp1, all.mergedPeaks = gr))
    }
}
.getReadLengthFromCigar <-function(cigar)
{
    if (is.na(cigar))
    {
        0
    }
    else if (length(grep("D", cigar) >0))
    {
        i <- substr(cigar, 1, nchar(cigar) - 1)
        sum(as.numeric(unlist(strsplit(gsub("^M", "", gsub("M|I|D|S", "M", i)),
            "M"))), na.rm = TRUE) - 2 * sum(as.numeric(lapply(unlist(strsplit(
            gsub("M|I|S", "M", i ), "M")), function(thisStr)
           {if (length(grep("D",thisStr)) >0)
               as.numeric(strsplit(thisStr, "D")[[1]][1]) else 0})))
    }      
    else
    { 
        i <- substr(cigar, 1, nchar(cigar) - 1)
        sum(as.numeric(unlist(strsplit(gsub("^M", "",
            gsub("M|I|D|S", "M", i)), "M"))), na.rm = TRUE)
    }
}
.getStrandedCoverage <-
function(gr, window.size = 20L, step = 10L,
   bg.window.size = 5000L, strand, min.reads = 10L)
{
    cvg <- coverage(gr)
    cvg <- Filter(length, cvg)
    observed <- runsum(cvg, k = window.size + 1, endrule = "constant")
    bg <- runsum(cvg, k = bg.window.size + 1, endrule = "constant")
    pos.value <- do.call(rbind, lapply(1:length(observed), function(i) {
        end <- cumsum(runLength(observed[[i]]))
        start <- end - runLength(observed[[i]]) + 1
        value <- runValue(observed[[i]])
        temp <- cbind(names(observed)[i], start, end, value)
        temp <- subset(temp, as.numeric(temp[,4]) >= min.reads)
        pos.bg <- as.numeric(temp[,2]) - ceiling((bg.window.size - window.size)/2)
        pos.bg[pos.bg < 1] <- 1
        bg.value <- as.data.frame(bg[[i]][pos.bg ])
        temp <- cbind(temp, bg.value)
        temp
    }))
    window.gr <- GRanges(IRanges(
        start = as.numeric(as.character(pos.value[,2])),
        end = as.numeric(as.character(pos.value[,3]))),
        seqnames = pos.value[,1], strand = Rle(strand, dim(pos.value)[1]),
        count = as.numeric(as.character(pos.value[,4])),
        bg = as.numeric(as.character(
            pos.value[,5])) / bg.window.size * window.size)
    window.gr
}
.locMaxPos <- function(data.ranges, window.size, step, min.reads)
{
    if (length(data.ranges) > 1)
    {
        total <- max(start(data.ranges))
        max.count <- numeric()
        i <- min(start(data.ranges))
        max.pos <- max.count
        j <- 0
        while (i < total){
            this.range <- data.ranges[start(data.ranges) >= i &
                start(data.ranges) < (i + window.size)]
            i <- i + step
            if (length(this.range) > 1)
            {
                max.ind <- this.range[which.max(as.numeric(
                    as.character(this.range$count)))]
                max.pos <- c(start(max.ind), max.pos)
                max.count <- c(max.ind$count, max.count)
            }
            else if (length(this.range) == 1)
            {
                max.pos <- c(start(this.range), max.pos)
                max.count <- c(as.numeric(as.character(this.range$count)),
                    max.count)
            }
            j <- min(which(start(data.ranges) > i))
            #print(paste("j = ", j))
            if (i < total & start(data.ranges)[j] >= (i + window.size))
            {
                i <- start(data.ranges)[j]
            }
            local.max.start <- unique(max.pos[max.count >= min.reads])
        }
    }
    else
    {
        local.max.start <- start(data.ranges[data.ranges$count >= min.reads])
    }
    #list(max.pos,  max.count)
    local.max.start 
}

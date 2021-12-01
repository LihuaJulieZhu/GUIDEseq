getUniqueCleavageEvents <-
    function(alignment.inputfile,
    umi.inputfile, 
    alignment.format = c("auto", "bam", "bed"),
    umi.header = FALSE,
    read.ID.col = 1,
    umi.col = 2,
    umi.sep = "\t",
    keep.chrM = FALSE,
    keep.R1only = TRUE,
    keep.R2only = TRUE, 
    concordant.strand = TRUE,
    max.paired.distance = 1000,
    min.mapping.quality = 30, 
    max.R1.len = 130,
    max.R2.len = 130,
    apply.both.max.len = FALSE,
    same.chromosome = TRUE, 
    distance.inter.chrom = -1,
    min.R1.mapped = 20, 
    min.R2.mapped = 20,
    apply.both.min.mapped = FALSE, 
    max.duplicate.distance = 0,
    umi.plus.R1start.unique = TRUE,
    umi.plus.R2start.unique = TRUE,
    n.cores.max = 6)
{ 
    if(!file.exists(alignment.inputfile))
        stop("alignment.inputfile is required, 
         please type ?getUniqueCleavageEvents for details!")
    if(!file.exists(umi.inputfile))
        stop("umi.input is required, please 
            type ?getUniqueCleavageEvents for details!")
### FIXME: could the UMI calculation be done in R, on the fly?
### FIXME: if not, it should be annotated in the read group of the BAM
    umi <- as.data.frame(fread(umi.inputfile, sep = umi.sep,
                               colClasses = "character", header = umi.header))
    alignment.format <- match.arg(alignment.format)
    if (alignment.format == "auto") {
        alignment.format <- tools::file_ext(alignment.inputfile)
    }
    if (alignment.format == "bed") {
        warning("BED alignment input is deprecated, please provide the BAM file")
        align <- importBEDAlignments(alignment.inputfile,
                                     min.mapping.quality, keep.chrM,
                                     keep.R1only, keep.R2only,
                                     min.R1.mapped, min.R2.mapped,
                                     apply.both.min.mapped,
                                     concordant.strand,
                                     same.chromosome,
                                     max.paired.distance,
                                     distance.inter.chrom,
                                     n.cores.max)
    } else {
        align <- importBAMAlignments(alignment.inputfile,
                                     min.mapping.quality, keep.chrM,
                                     keep.R1only, keep.R2only,
                                     min.R1.mapped, min.R2.mapped,
                                     apply.both.min.mapped,
                                     concordant.strand,
                                     same.chromosome,
                                     max.paired.distance,
                                     distance.inter.chrom
                                     )
    }
    pdf(paste0(gsub(".sorted", "", gsub(".bam", "", basename(alignment.inputfile))),
        "AlignmentWidthDistribution.pdf"))
    par(mfrow = c(1,2))
    hist(align$qwidth.first, xlab = "R1 aligned read width", main = "")
    hist(align$qwidth.last, xlab = "R2 aligned read width", main = "")
    dev.off() 
    if (length(umi) < read.ID.col || length(umi) < umi.col)
    {
        stop("umi input file must contain at least two columns, 
            one is the header of the read, specified by read.ID.col,
            the other column is the umi sequence column, specified by umi.col")
    } 
    else
    {
        umi <- umi[, c(read.ID.col, umi.col)]
        colnames(umi) <- c("readName", "UMI")
        umi$readName <- gsub("^@", "", umi$readName)
        if (apply.both.max.len)
            align <- subset(align, qwidth.first <= max.R1.len &
                qwidth.last <= max.R2.len)
        else
            align <- subset(align,
                            (qwidth.first > 0L & qwidth.first <= max.R1.len) |
                                (qwidth.last > 0L & qwidth.last <= max.R2.len))
        align.umi <- merge(align, umi)

### plus means R2 on plus strand
        R2.good.len <- subset(align.umi, qwidth.last <= max.R2.len & 
                                  qwidth.last >= min.R2.mapped)
        R2.umi.plus <- subset(R2.good.len, strand.last == "+")
        R2.umi.minus <- subset(R2.good.len, strand.last == "-")

        R1.good.len <- subset(align.umi, qwidth.first <= max.R1.len & 
                                  qwidth.first >= min.R1.mapped)
        R1.umi.plus <- subset(R1.good.len, strand.first == "-" &
                                  !readName %in% R2.umi.plus$readName)
        R1.umi.minus <- subset(R1.good.len, strand.first == "+" &
                                   !readName %in% R2.umi.minus$readName)
        
        unique.umi.plus.R2 <-
            unique(R2.umi.plus[, c("seqnames.last", "seqnames.first", 
                                   "strand.last", "strand.first", "start.last",
                                   "end.first", "UMI")])
        unique.umi.minus.R2 <-
            unique(R2.umi.minus[, c("seqnames.last", "seqnames.first",  
                                    "strand.last", "strand.first",
                                    "end.last", "start.first", "UMI")])
        unique.umi.plus.R1 <-
            unique(R1.umi.plus[, c("seqnames.last", "seqnames.first", 
                                   "strand.last", "strand.first",
                                   "start.first", "start.last", "UMI")])
        unique.umi.minus.R1 <-
            unique(R1.umi.minus[, c("seqnames.last", "seqnames.first", 
                                    "strand.last", "strand.first",
                                    "end.first", "end.last", "UMI")])
        plus.cleavage.R2 <-
            unique.umi.plus.R2[, c("seqnames.last", "start.last")]
        plus.cleavage.R1 <-
            unique.umi.plus.R1[, c("seqnames.first", "start.first")]
        minus.cleavage.R2 <-
            unique.umi.minus.R2[, c("seqnames.last", "end.last")]
        minus.cleavage.R1 <-
            unique.umi.minus.R1[, c("seqnames.first", "end.first")]
        colnames(plus.cleavage.R1) <- c("seqnames", "start")
        colnames(plus.cleavage.R2) <- c("seqnames", "start")
        colnames(minus.cleavage.R1) <- c("seqnames", "start")
        colnames(minus.cleavage.R2) <- c("seqnames", "start")
        plus.cleavage <- rbind(plus.cleavage.R2, plus.cleavage.R1)
        minus.cleavage <- rbind(minus.cleavage.R2, minus.cleavage.R1)
        plus.cleavage <- cbind(plus.cleavage, strand = "+")
        minus.cleavage <- cbind(minus.cleavage, strand = "-")
        colnames(plus.cleavage) <- c("seqnames", "start", "strand")
        colnames(minus.cleavage) <- c("seqnames", "start", "strand")
        unique.umi.both <- rbind(plus.cleavage, minus.cleavage)

        R1.umi.plus <- R1.umi.plus[, c("seqnames.first",
                                   "strand.first",
                                   "start.first", "UMI")]
        R1.umi.minus <- R1.umi.minus[, c("seqnames.first",
                                    "strand.first",
                                    "end.first", "UMI")]
        R2.umi.plus <- R2.umi.plus[, c("seqnames.last",
                                    "strand.last",
                                    "start.last", "UMI")] 
        R2.umi.minus <- R2.umi.minus[, c("seqnames.last",
                                    "strand.last", 
                                    "end.last", "UMI")]

        colnames(R1.umi.plus) <- c("seqnames", "strand", "start", "UMI")
        colnames(R1.umi.minus) <- c("seqnames", "strand", "start", "UMI")
        colnames(R2.umi.plus) <- c("seqnames", "strand", "start", "UMI")
        colnames(R2.umi.minus) <- c("seqnames", "strand", "start", "UMI")

	R1.umi.plus.summary <- unique(add_count(R1.umi.plus, seqnames, strand, start, UMI))
        R1.umi.minus.summary <- unique(add_count(R1.umi.minus, seqnames, strand, start, UMI))
        R2.umi.plus.summary <- unique(add_count(R2.umi.plus, seqnames, strand, start, UMI))
        R2.umi.minus.summary <- unique(add_count(R2.umi.minus, seqnames, strand, start, UMI))
         
        list(cleavage.gr = GRanges(IRanges(start=unique.umi.both[,2], width=1),
            seqnames=unique.umi.both[,1], strand = unique.umi.both[,3], 
            total=rep(1, dim(unique.umi.both)[1])), 
            unique.umi.plus.R2 = unique.umi.plus.R2, 
            unique.umi.minus.R2 = unique.umi.minus.R2,
            unique.umi.plus.R1 = unique.umi.plus.R1, 
            unique.umi.minus.R1 = unique.umi.minus.R1, 
            align.umi = align.umi,
            umi.count.summary = rbind(R1.umi.plus.summary, 
              R1.umi.minus.summary,
              R2.umi.plus.summary,
              R2.umi.minus.summary)
	) 
    }
}

rbind_dodge <- function(x, y,
                        xSuffix = deparse(substitute(x)),
                        ySuffix = deparse(substitute(y)),
                        common = character())
{
    naDF <- function(cols) {
        as.data.frame(setNames(as.list(rep(NA, length(cols))), cols))
    }

    dodgeDF <- function(df, suffix) {
        only <- setdiff(colnames(df), common)
        dodged <- df[only]
        colnames(dodged) <- paste0(only, suffix)
        dodged
    }

    xDodged <- dodgeDF(x, xSuffix)    
    yDodged <- dodgeDF(y, ySuffix)
    
    rbind(cbind(x[common], xDodged, naDF(colnames(yDodged))),
          cbind(y[common], naDF(colnames(xDodged)), yDodged))
}


importBAMAlignments <- function(file,
                                min.mapping.quality = 30L,
                                keep.chrM = FALSE,
                                keep.R1only = TRUE,
                                keep.R2only = TRUE, 
                                min.R1.mapped = 30L, 
                                min.R2.mapped = 30L,
                                apply.both.min.mapped = FALSE,
                                concordant.strand = TRUE,
                                same.chromosome = TRUE,
                                max.paired.distance = 1000L,
                                distance.inter.chrom = -1L
                                )
{
    param <- ScanBamParam(mapqFilter=min.mapping.quality, what=c("flag", "mapq"))
    gal <- readGAlignmentsList(BamFile(file, asMates=TRUE), param=param,
                               use.names=TRUE)
    my.pairs <- tryCatch(( as(gal, "GAlignmentPairs")),
          error=function(e) { 
             readGAlignmentPairs(BamFile(file), use.names = TRUE, 
                param = param, strandMode = 1)})

    if (!keep.chrM)
        my.pairs <- my.pairs[!seqnames(my.pairs) %in% c("chrM", "chrMT", "M"),]
    if (apply.both.min.mapped) {
        my.pairs <- my.pairs[width(first(my.pairs)) >= min.R1.mapped &
                       width(last(my.pairs)) >= min.R2.mapped]
    } else {
        my.pairs <- my.pairs[width(first(my.pairs)) >= min.R1.mapped |
                       width(last(my.pairs)) >= min.R2.mapped]
    }

    if (concordant.strand) {
        my.pairs <- my.pairs[strand(my.pairs) != "*"]
    }

    if (same.chromosome) {
        my.pairs <- my.pairs[!is.na(seqnames(my.pairs))]
    }
    
    distance <- ifelse(is.na(seqnames(my.pairs)), distance.inter.chrom,
                       ## does not yield negative
                       ## width(pgap(ranges(first(my.pairs)), ranges(last(my.pairs))))
                       ifelse(strand(my.pairs) == "+",
                              (start(last(my.pairs)) - end(first(my.pairs))), 
                              (start(first(my.pairs)) - end(last(my.pairs)))))
    mcols(my.pairs)$distance <- distance

    unpaired <- unlist(gal[mcols(gal)$mate_status == "unmated"])

    mcols(unpaired)$readName <- names(unpaired)
    names(unpaired) <- NULL
    
    first <- bamFlagTest(mcols(unpaired)$flag, "isFirstMateRead")
    firstUnpaired <- unpaired[first & keep.R1only]
    firstUnpaired <- firstUnpaired[width(firstUnpaired) >= min.R1.mapped]
    lastUnpaired <- unpaired[!first & keep.R2only]
    lastUnpaired <- lastUnpaired[width(lastUnpaired) >= min.R2.mapped]

    mcols(unpaired)$flag <- NULL
    unpairedDF <- rbind_dodge(as.data.frame(firstUnpaired),
                              as.data.frame(lastUnpaired),
                              ".first", ".last", "readName")
    unpairedDF$distance <- NA_integer_
    
    pairedDF <- as.data.frame(my.pairs)
    pairedDF$readName <- rownames(pairedDF)
    rownames(pairedDF) <- NULL
    
    df <- rbind(pairedDF[colnames(unpairedDF)], unpairedDF)
    df$njunc.first <- df$njunc.last <- NULL
    
    df
}

importBEDAlignments <- function(file,
                                min.mapping.quality = 30L,
                                keep.chrM = FALSE,
                                keep.R1only = TRUE,
                                keep.R2only = TRUE, 
                                min.R1.mapped = 30L, 
                                min.R2.mapped = 30L,
                                apply.both.min.mapped = FALSE,
                                concordant.strand = TRUE,
                                same.chromosome = TRUE,
                                max.paired.distance = 1000L,
                                distance.inter.chrom = -1L,
                                n.cores = 6)
{
    align <- as.data.frame(fread(file, sep = "\t",
                                  header = FALSE, 
                                 colClasses =
                                     c("character", "integer",
                                       "integer", "character", "integer",
                                       "character", "character")))
    colnames(align) <- c("seqnames", "start", "end", "name", "mapq",
                         "strand", "cigar")
    align <- align[align$mapq >= min.mapping.quality,]
    if (!keep.chrM)
        align <- align[!align$seqnames %in% c("chrM", "chrMT", "M"),]
    sidePos <- nchar(align$name)
    readSide <- as.integer(substring(align$name, sidePos, sidePos))
    align$readName <- substring(align$name, 1L, sidePos - 2L)
    align$name <- NULL
    R1 <- align[readSide == 1L, ]
    R2 <- align[readSide == 2L, ]
    all <- merge(R1, R2, by="readName", all.x = keep.R1only, 
                 all.y = keep.R2only, suffixes = c(".first", ".last"))
    if (apply.both.min.mapped)
        all <- subset(all, (is.na(all$end.first) | 
                           (all$end.first - all$start.first) >= min.R1.mapped) &
                             (is.na(all$end.last) | 
                             (all$end.last - all$start.last) >= min.R2.mapped))
    else
        all <- subset(all, (all$end.first - all$start.first) >= min.R1.mapped |
                          (all$end.last - all$start.last) >= min.R2.mapped )
    if (concordant.strand)
        all <- subset(all, is.na(all$strand.first) | is.na(all$strand.last) |
                          all$strand.last != all$strand.first)
    all$start.first <- all$start.first + 1L # BED is 0-based
    all$start.last <- all$start.last + 1L
    if (same.chromosome)
        all <- subset(all, is.na(all$seqnames.first) | is.na(all$seqnames.last) |
                          all$seqnames.first == all$seqnames.last)
    if (keep.R1only && !keep.R2only)
    {
        all <- subset(all, !is.na(all$seqnames.last))
    }
    else if (keep.R2only && !keep.R1only)
    {
        all <- subset(all, !is.na(all$seqnames.first))
    }
    else if (!keep.R1only && !keep.R2only)
    {
        all <- subset(all, !is.na(all$seqnames.first) &
                          !is.na(all$seqnames.last))
    }
    distance <- ifelse(all$strand.last == "-", (all$start.last - all$end.first), 
                       (all$start.first - all$end.last))
    distance[!is.na(all$seqnames.first) & !is.na(all$seqnames.last) &
                 all$seqnames.first != all$seqnames.last] <- distance.inter.chrom
    all <- cbind(all, distance)
    all <- subset(all, is.na(distance) | distance <=  max.paired.distance)
    n.cores <- min(n.cores, detectCores() -1 )
    unique.cigar <- unique(c(all$cigar.first, all$cigar.last))
    if (n.cores > 1)
    {
        cl <- makeCluster(n.cores)
        unique.base.kept <- cbind(cigar = unique.cigar,
                              base.kept = unlist(parLapply(cl, unique.cigar, 
                              .getReadLengthFromCigar)))
        stopCluster(cl)
    }
    else
    {
        unique.base.kept <- cbind(cigar = unique.cigar,
                              base.kept = unlist(lapply(unique.cigar,
                              .getReadLengthFromCigar)))
    }
    qwidth.first <- as.numeric(
        unique.base.kept[match(all$cigar.first, unique.base.kept[,1]),2])
    qwidth.last <- as.numeric(
        unique.base.kept[match(all$cigar.last, unique.base.kept[,1]),2])
    cbind(all, qwidth.first, qwidth.last)
}

getUniqueCleavageEvents <-
    function(alignment.inputfile,
    umi.inputfile, 
    alignment.format = "bed",
    umi.header = FALSE,
    read.ID.col = 1,
    umi.col = 2,
    umi.sep = "\t",
    keep.R1only = TRUE,
    keep.R2only = TRUE, 
    paired.direction = "opposite.strand",
    max.paired.distance = 1000,
    min.mapping.quality = 30, 
    max.R1.len = 130,
    max.R2.len = 130,
    apply.both.max.len = FALSE,
    same.chromosome = TRUE, 
    distance.inter.chrom = -1,
    min.R1.mapped = 30, 
    min.R2.mapped = 30,
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
    umi <- as.data.frame(fread(umi.inputfile, sep = umi.sep,
        colClasses = "character", header = umi.header))
    align <- as.data.frame(fread(alignment.inputfile, sep = "\t",
        header = FALSE, 
        colClasses = c("character", "numeric", "numeric", 
        "character", "numeric", "character", "character")))
    if (dim(umi)[1] < read.ID.col || dim(umi)[1] < umi.col )
    {
        stop("umi input file must contain at least two columns, 
            one is the header of the read, specified by read.ID.col,
            the other column is the umi sequence column, specified by umi.col")
    } 
    else
    {
        umi <- umi[, c(read.ID.col, umi.col)]
        colnames(umi) <- c("readName", "UMI")
        align <- subset(align, align[,5] >= min.mapping.quality)
        n.cores <- detectCores()
        n.cores <- min(n.cores, n.cores.max)
        if (n.cores > 1)
        {
            cl <- makeCluster(n.cores)
            reads <- do.call(rbind, parLapply(cl, align[,4], function(i) {
                unlist(strsplit(i, "/")) }))
            stopCluster(cl)
        }
        else
        {
            reads <- do.call(rbind, lapply(align[,4], function(i) {
                unlist(strsplit(i, "/")) }))
        }
        align <- cbind(align, reads)
        align <- align[, -4]
        colnames(align) <- c("chr", "start", "end","mapping.qual",
            "strand", "cigar","readName", "readSide")
        R1 <- align[align[,8] == 1, ]
        R2 <- align[align[,8] == 2, ]
        all <- merge(R1, R2, by="readName", all.x = keep.R1only, 
            all.y = keep.R2only)
        if (apply.both.min.mapped)
            all <- subset(all, (is.na(all$end.x) | 
                (all$end.x - all$start.x) >= min.R1.mapped) &
                (is.na(all$end.y) | 
                (all$end.y - all$start.y) >= min.R2.mapped))
        else
            all <- subset(all, (all$end.x - all$start.x) >= min.R1.mapped |
                (all$end.y - all$start.y) >= min.R2.mapped )
        if (paired.direction == "opposite.strand")
            all <- subset(all, is.na(all$strand.x) | is.na(all$strand.y) |
                all$strand.y != all$strand.x)
        if (same.chromosome)
            all <- subset(all, is.na(all$chr.x) | is.na(all$chr.y) |
                all$chr.x == all$chr.y)
        if (keep.R1only && !keep.R2only)
        {
            all <- subset(all, !is.na(all$chr.y))
        }
        else if (keep.R2only && !keep.R1only)
        {
            all <- subset(all, !is.na(all$chr.x))
        }
        else if (!keep.R1only && !keep.R2only)
        {
            all <- subset(all, !is.na(all$chr.x) & !is.na(all$chr.y))
        }
        distance <- ifelse(all$strand.y == "-", (all$start.y - all$end.x), 
            (all$start.x - all$end.y))
        distance[!is.na(all$chr.x) & !is.na(all$chr.y) &
            all$chr.x != all$chr.y] <- distance.inter.chrom
        all <- cbind(all, distance)
        all <- subset(all, is.na(distance) | distance <=  max.paired.distance)
        all[,1] <- gsub("@", "", all[,1])
        umi[,1] <- gsub("@", "", umi[,1])
        unique.cigar <- unique(c(all$cigar.x, all$cigar.y))
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
        R1.base.kept <- as.numeric(
            unique.base.kept[match(all$cigar.x, unique.base.kept[,1]),2])
        R2.base.kept <- as.numeric(
            unique.base.kept[match(all$cigar.y, unique.base.kept[,1]),2])
        all <- cbind(all, R1.base.kept, R2.base.kept)
        if (apply.both.max.len)
            all <- subset(all, R1.base.kept <= max.R1.len &
                R2.base.kept <= max.R2.len)
        else
            all <- subset(all, (R1.base.kept > 0 & R1.base.kept <= max.R1.len) |
                (R2.base.kept >0 & R2.base.kept <= max.R2.len))
        all.umi <- merge(all, umi)
        ### plus means R2 on plus strand
        R2.umi.plus <- subset(all.umi, !is.na(all.umi$strand.y) & 
            all.umi$strand.y == "+" &
            all.umi$R2.base.kept <= max.R2.len & 
            all.umi$R2.base.kept >= min.R2.mapped)
        R2.umi.minus <- subset(all.umi, !is.na(all.umi$strand.y) & 
            all.umi$strand.y == "-" &
            all.umi$R2.base.kept <= max.R2.len & 
            all.umi$R2.base.kept >= min.R2.mapped)
        R1.umi.plus <- subset(all.umi, !is.na(all.umi$strand.x) & 
            all.umi$strand.x == "-" &
            all.umi$R1.base.kept <= max.R1.len & 
            all.umi$R1.base.kept >= min.R1.mapped &
            !all.umi$readName %in% R2.umi.plus$readName)
        R1.umi.minus <- subset(all.umi, !is.na(all.umi$strand.x) & 
            all.umi$strand.x == "+" &
            all.umi$R1.base.kept <= max.R1.len & 
            all.umi$R1.base.kept >= min.R1.mapped &
            !all.umi$readName %in% R2.umi.minus$readName) 
        unique.umi.plus.R2 <- unique(R2.umi.plus[, c("chr.y", "chr.x", 
            "strand.y", "strand.x", "start.y", "end.x", "UMI")])
        unique.umi.minus.R2 <- unique(R2.umi.minus[, c("chr.y", "chr.x",  
            "strand.y", "strand.x", "end.y", "start.x", "UMI")])
        unique.umi.plus.R1 <- unique(R1.umi.plus[, c("chr.y", "chr.x", 
            "strand.y", "strand.x", "start.x", "start.y", "UMI")])
        unique.umi.minus.R1 <- unique(R1.umi.minus[, c("chr.y", "chr.x", 
            "strand.y", "strand.x", "end.x", "end.y", "UMI")])
        plus.cleavage.R2 <- unique.umi.plus.R2[, c("chr.y", "start.y")]
        plus.cleavage.R1 <- unique.umi.plus.R1[, c("chr.x", "start.x")]
        minus.cleavage.R2 <- unique.umi.minus.R2[, c("chr.y", "end.y")]
        minus.cleavage.R1 <- unique.umi.minus.R1[, c("chr.x", "end.x")]
        colnames(plus.cleavage.R1) <- c("chr", "start")
        colnames(plus.cleavage.R2) <- c("chr", "start")
        colnames(minus.cleavage.R1) <- c("chr", "start")
        colnames(minus.cleavage.R2) <- c("chr", "start")
        plus.cleavage <- rbind(plus.cleavage.R2, plus.cleavage.R1)
        minus.cleavage <- rbind(minus.cleavage.R2, minus.cleavage.R1)
        plus.cleavage <- cbind(plus.cleavage, strand = "+")
        minus.cleavage <- cbind(minus.cleavage, strand = "-")
        colnames(plus.cleavage) <- c("chr", "start", "strand")
        colnames(minus.cleavage) <- c("chr", "start", "strand")
        unique.umi.both <- rbind(plus.cleavage, minus.cleavage)
        list(cleavage.gr = GRanges(IRanges(start=unique.umi.both[,2], width=1),
            seqnames=unique.umi.both[,1], strand = unique.umi.both[,3], 
            total=rep(1, dim(unique.umi.both)[1])), 
            unique.umi.plus.R2 = unique.umi.plus.R2, 
            unique.umi.minus.R2 = unique.umi.minus.R2,
            unique.umi.plus.R1 = unique.umi.plus.R1, 
            unique.umi.minus.R1 = unique.umi.minus.R1, 
            all.umi = all.umi) 
    }
}

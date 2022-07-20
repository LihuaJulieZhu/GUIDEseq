#' Using UMI sequence to obtain the starting sequence library
#'
#' PCR amplification often leads to biased representation of the starting
#' sequence population. To track the sequence tags present in the initial
#' sequence library, a unique molecular identifier (UMI) is added to the 5
#' prime of each sequence in the staring library. This function uses the UMI
#' sequence plus the first few sequence from R1 reads to obtain the starting
#' sequence library.
#'
#'
#' @param alignment.inputfile The alignment file. Currently supports bed output
#' file with CIGAR information. Suggest run the workflow binReads.sh, which
#' sequentially runs barcode binning, adaptor removal, alignment to genome,
#' alignment quality filtering, and bed file conversion. Please download the
#' workflow function and its helper scripts at
#' http://mccb.umassmed.edu/GUIDE-seq/binReads/
#' @param umi.inputfile A text file containing at least two columns, one is the
#' read identifier and the other is the UMI or UMI plus the first few bases of
#' R1 reads. Suggest use getUMI.sh to generate this file. Please download the
#' script and its helper scripts at http://mccb.umassmed.edu/GUIDE-seq/getUMI/
#' @param alignment.format The format of the alignment input file. Currently
#' only bam and bed file format is supported. BED format will be deprecated
#' soon.
#' @param umi.header Indicates whether the umi input file contains a header
#' line or not.  Default to FALSE
#' @param read.ID.col The index of the column containing the read identifier in
#' the umi input file, default to 1
#' @param umi.col The index of the column containing the umi or umi plus the
#' first few bases of sequence from the R1 reads, default to 2
#' @param umi.sep column separator in the umi input file, default to tab
#' @param keep.chrM Specify whether to include alignment from chrM. Default
#' FALSE
#' @param keep.R1only Specify whether to include alignment with only R1 without
#' paired R2.  Default TRUE
#' @param keep.R2only Specify whether to include alignment with only R2 without
#' paired R1.  Default TRUE
#' @param concordant.strand Specify whether the R1 and R2 should be aligned to
#' the same strand or opposite strand. Default opposite.strand (TRUE)
#' @param max.paired.distance Specify the maximum distance allowed between
#' paired R1 and R2 reads.  Default 1000 bp
#' @param min.mapping.quality Specify min.mapping.quality of acceptable
#' alignments
#' @param max.R1.len The maximum retained R1 length to be considered for
#' downstream analysis, default 130 bp. Please note that default of 130 works
#' well when the read length 150 bp. Please also note that retained R1 length
#' is not necessarily equal to the mapped R1 length
#' @param max.R2.len The maximum retained R2 length to be considered for
#' downstream analysis, default 130 bp. Please note that default of 130 works
#' well when the read length 150 bp. Please also note that retained R2 length
#' is not necessarily equal to the mapped R2 length
#' @param apply.both.max.len Specify whether to apply maximum length
#' requirement to both R1 and R2 reads, default FALSE
#' @param same.chromosome Specify whether the paired reads are required to
#' align to the same chromosome, default TRUE
#' @param distance.inter.chrom Specify the distance value to assign to the
#' paired reads that are aligned to different chromosome, default -1
#' @param min.R1.mapped The maximum mapped R1 length to be considered for
#' downstream analysis, default 30 bp.
#' @param min.R2.mapped The maximum mapped R2 length to be considered for
#' downstream analysis, default 30 bp.
#' @param apply.both.min.mapped Specify whether to apply minimum mapped length
#' requirement to both R1 and R2 reads, default FALSE
#' @param max.duplicate.distance Specify the maximum distance apart for two
#' reads to be considered as duplicates, default 0. Currently only 0 is
#' supported
#' @param umi.plus.R1start.unique To specify whether two mapped reads are
#' considered as unique if both containing the same UMI and same alignment
#' start for R1 read, default TRUE.
#' @param umi.plus.R2start.unique To specify whether two mapped reads are
#' considered as unique if both containing the same UMI and same alignment
#' start for R2 read, default TRUE.
#' @param min.umi.count To specify the minimum count for a umi to be included
#' in the subsequent analysis.  Please adjust it to a higher number for deeply
#' sequenced library and vice versa.
#' @param max.umi.count To specify the maximum count for a umi to be included
#' in the subsequent analysis.  Please adjust it to a higher number for deeply
#' sequenced library and vice versa.
#' @param min.read.coverage To specify the minimum coverage for a read UMI
#' combination to be included in the subsequent analysis.  Please note that
#' this is different from min.umi.count which is less stringent.
#' @param n.cores.max Indicating maximum number of cores to use in multi core
#' mode, i.e., parallel processing, default 6. Please set it to 1 to disable
#' multicore processing for small dataset.
#' @param outputDir output Directory to save the figures
#' @return \item{cleavage.gr }{Cleavage sites with one site per UMI as GRanges
#' with metadata column total set to 1 for each range}
#' \item{unique.umi.plus.R2}{a data frame containing unique cleavage site from
#' R2 reads mapped to plus strand with the following columns seqnames
#' (chromosome) start (cleavage site) strand UMI (unique molecular identifier
#' (umi) or umi with the first few bases of R1 read) UMI read duplication level
#' (min.read.coverage can be used to remove UMI-read with very low coverage) }
#' \item{unique.umi.minus.R2}{a data frame containing unique cleavage site from
#' R2 reads mapped to minus strand with the same columns as unique.umi.plus.R2
#' } \item{unique.umi.plus.R1}{a data frame containing unique cleavage site
#' from R1 reads mapped to minus strand without corresponding R2 reads mapped
#' to the plus strand, with the same columns as unique.umi.plus.R2 }
#' \item{unique.umi.minus.R1}{a data frame containing unique cleavage site from
#' R1 reads mapped to plus strand without corresponding R2 reads mapped to the
#' minus strand, with the same columns as unique.umi.plus.R2 } \item{all.umi}{a
#' data frame containing all the mapped reads with the following columns.
#' readName (read ID), chr.x (chromosome of readSide.x/R1 read), start.x (start
#' of eadSide.x/R1 read), end.x (end of eadSide.x/R1 read), mapping.qual.x
#' (mapping quality of readSide.x/R1 read), strand.x (strand of readSide.x/R1
#' read), cigar.x (CIGAR of readSide.x/R1 read) , readSide.x (1/R1), chr.y
#' (chromosome of readSide.y/R2 read) start.y (start of readSide.y/R2 read),
#' end.y (end of readSide.y/R2 read), mapping.qual.y (mapping quality of
#' readSide.y/R2 read), strand.y (strand of readSide.y/R2 read), cigar.y (CIGAR
#' of readSide.y/R2 read), readSide.y (2/R2) R1.base.kept (retained R1 length),
#' R2.base.kept (retained R2 length), distance (distance between mapped R1 and
#' R2), UMI (unique molecular identifier (umi) or umi with the first few bases
#' of R1 read) }
#' @author Lihua Julie Zhu
#' @seealso getPeaks
#' @references Shengdar Q Tsai and J Keith Joung et al. GUIDE-seq enables
#' genome-wide profiling of off-target cleavage by CRISPR-Cas nucleases. Nature
#' Biotechnology 33, 187 to 197 (2015)
#' @keywords misc
#' @examples
#'
#'     if(interactive())
#'     {
#'         umiFile <- system.file("extdata", "UMI-HEK293_site4_chr13.txt",
#'            package = "GUIDEseq")
#'         alignFile <- system.file("extdata","bowtie2.HEK293_site4_chr13.sort.bam" ,
#'             package = "GUIDEseq")
#'         cleavages <- getUniqueCleavageEvents(
#'             alignment.inputfile = alignFile , umi.inputfile = umiFile,
#'             n.cores.max = 1)
#'         names(cleavages)
#'         #output a summary of duplicate counts for sequencing saturation assessment
#'         table(cleavages$umi.count.summary$n)
#'     }
#' @importFrom stats setNames
#' @importFrom methods as
#' @rawNamespace import(S4Vectors, except=c(fold, values, rename))
#' @rawNamespace import(IRanges, except=values)
#' @rawNamespace import(BSgenome, except=export)
#' @rawNamespace import(GenomicRanges, except=values)
#' @rawNamespace import(BiocGenerics, except=c(var, sd))
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics hist par
#' @importFrom GenomicAlignments readGAlignmentsList first
#' last readGAlignmentPairs
#' @importFrom data.table fread
#' @importFrom parallel makeCluster stopCluster detectCores
#' parLapply
#' @importFrom Rsamtools ScanBamParam BamFile bamFlagTest
#' @importFrom tools file_ext
#' @importFrom dplyr select mutate add_count filter '%>%'

#' @export getUniqueCleavageEvents
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
    max.duplicate.distance = 0L,
    umi.plus.R1start.unique = TRUE,
    umi.plus.R2start.unique = TRUE,
    min.umi.count = 5L,
    max.umi.count = 100000L,
    min.read.coverage = 1L,
    n.cores.max = 6,
    outputDir)
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
        alignment.format <- file_ext(alignment.inputfile)
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
    if (missing(outputDir))
    {
        outputDir <- getwd()
    }
    pdf(file.path(outputDir, paste0(gsub(".sorted", "", gsub(".bam", "",
                    basename(alignment.inputfile))),
        "AlignmentWidthDistribution.pdf")))
    par(mfrow = c(1,2))
    hist(align$width.first, xlab = "R1 aligned read width", main = "")
    hist(align$width.last, xlab = "R2 aligned read width", main = "")
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
            align <- subset(align, width.first <= max.R1.len &
                width.last <= max.R2.len)
        else
            align <- subset(align,
                            (width.first > 0L & width.first <= max.R1.len) |
                                (width.last > 0L & width.last <= max.R2.len))
        align.umi <- merge(align, umi)
        all.ind <- seq(dim(align.umi)[1])
        align.umi <- align.umi[setdiff(all.ind, grep("N", align.umi$UMI)), ]
        temp <- as.data.frame(table(align.umi$UMI))
	align.umi <- align.umi[align.umi$UMI %in% temp[temp[,2] >= min.umi.count &
                   temp[,2] <= max.umi.count,1],]
        pdf(file.path(outputDir,paste0(gsub(".sorted", "", gsub(".bam", "",
                   basename(alignment.inputfile))),
              "AlignmentWidthDistributionUmiWithoutN.pdf")))
              par(mfrow = c(1,2))
              hist(align.umi$width.first, xlab = "R1 aligned read width", main = "")
              hist(align.umi$width.last, xlab = "R2 aligned read width", main = "")
        dev.off()
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

        unique.umi.plus.R2 <- R2.umi.plus %>%
            select(seqnames.last, seqnames.first,
                     strand.last, strand.first,
                     start.last, end.first, UMI) %>%
            add_count(seqnames.last, seqnames.first,
                     strand.last, strand.first,
                     start.last, end.first, UMI) %>%
            unique %>%
            filter(n >= min.read.coverage)

        unique.umi.minus.R2 <- R2.umi.minus %>%
            select(seqnames.last, seqnames.first,
                     strand.last, strand.first,
                     end.last, start.first, UMI) %>%
            add_count(seqnames.last, seqnames.first,
                     strand.last, strand.first,
                     end.last, start.first, UMI) %>%
            unique %>%
            filter(n >= min.read.coverage)

        unique.umi.plus.R1 <- R1.umi.plus %>%
            select(seqnames.last, seqnames.first,
                     strand.last, strand.first,
                     start.first, start.last, UMI) %>%
            add_count(seqnames.last, seqnames.first,
                     strand.last, strand.first,
                     start.first, start.last, UMI) %>%
            unique %>%
            filter(n >= min.read.coverage)

        unique.umi.minus.R1 <- R1.umi.minus  %>%
            select(seqnames.last, seqnames.first,
                     strand.last, strand.first,
                     end.first, end.last, UMI) %>%
            add_count(seqnames.last, seqnames.first,
                     strand.last, strand.first,
                     end.first, end.last, UMI) %>%
            unique %>%
            filter(n >= min.read.coverage)

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

        res <- list(cleavage.gr = GRanges(IRanges(start=unique.umi.both[,2], width=1),
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
        #saveRDS(res, file = paste0(gsub(".sorted", "",
        #      gsub(".bam", "", basename(alignment.inputfile))),
        #      "CleavageSitesWithUMI.RDS"))
        res
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

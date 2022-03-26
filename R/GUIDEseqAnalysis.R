GUIDEseqAnalysis <- function(alignment.inputfile,
    umi.inputfile,
    alignment.format = c("auto", "bam", "bed"),
    umi.header = FALSE,
    read.ID.col = 1L,
    umi.col = 2L,
    umi.sep = "\t",
    BSgenomeName,
    gRNA.file,
    outputDir,
    n.cores.max = 1L,
    keep.chrM = FALSE,
    keep.R1only = TRUE,
    keep.R2only = TRUE,
    concordant.strand = TRUE,
    max.paired.distance = 1000L,
    min.mapping.quality = 30L,
    max.R1.len = 130L,
    max.R2.len = 130L,
    min.umi.count = 5L,
    max.umi.count = 100000L,
    min.read.coverage = 1L,
    apply.both.max.len = FALSE,
    same.chromosome = TRUE,
    distance.inter.chrom = -1L,
    min.R1.mapped = 20L,
    min.R2.mapped = 20L,
    apply.both.min.mapped = FALSE,
    max.duplicate.distance = 0L,
    umi.plus.R1start.unique = TRUE,
    umi.plus.R2start.unique = TRUE,
    window.size = 20L,
    step = 20L,
    bg.window.size = 5000L,
    min.reads = 5L,
    min.reads.per.lib = 1L,
    min.peak.score.1strandOnly = 5L,
    min.SNratio = 2,
    maxP = 0.01,
    stats = c("poisson", "nbinom"),
    p.adjust.methods =
        c( "none", "BH", "holm", "hochberg", "hommel", "bonferroni",
        "BY", "fdr"),
    distance.threshold = 40L,
    max.overlap.plusSig.minusSig = 30L,
    plus.strand.start.gt.minus.strand.end = TRUE,
    keepPeaksInBothStrandsOnly = TRUE,
    gRNA.format = "fasta",
    overlap.gRNA.positions = c(17,18),
    upstream = 20L,
    downstream = 20L,
    PAM.size = 3L,
    gRNA.size = 20L,
    PAM = "NGG",
    PAM.pattern = "NNN$",
    max.mismatch = 6L,
    allowed.mismatch.PAM = 2L,
    overwrite = TRUE,
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,0.804, 0.685, 0.583),
    orderOfftargetsBy = c("peak_score", "predicted_cleavage_score", "n.mismatch"),
    descending = TRUE,
    keepTopOfftargetsOnly = TRUE,
    keepTopOfftargetsBy = c("predicted_cleavage_score", "n.mismatch"),
    scoring.method = c("Hsu-Zhang", "CFDscore"),
        subPAM.activity = hash( AA =0,
          AC =   0,
          AG = 0.259259259,
          AT = 0,
          CA = 0,
          CC = 0,
          CG = 0.107142857,
          CT = 0,
          GA = 0.069444444,
          GC = 0.022222222,
          GG = 1,
          GT = 0.016129032,
          TA = 0,
          TC = 0,
          TG = 0.038961039,
          TT = 0),
     subPAM.position = c(22, 23),
     PAM.location = "3prime",
     mismatch.activity.file = system.file("extdata",
         "NatureBiot2016SuppTable19DoenchRoot.csv",
         package = "CRISPRseek"),
     txdb, orgAnn
)
{
    alignment.format <- match.arg(alignment.format)
    orderOfftargetsBy <- match.arg(orderOfftargetsBy)
    keepTopOfftargetsBy <- match.arg(keepTopOfftargetsBy)
    message("Remove duplicate reads ...\n")
    if (missing(BSgenomeName))
    {
        stop("genome is required parameter,
            please pass in either a BSgenome object!")
    }
    if (!class(BSgenomeName)== "BSgenome")
    {
        stop("BSgenomeName must be a BSgenome object!")
    }
    if (missing(alignment.inputfile) || missing(umi.inputfile) ||
        missing(gRNA.file))
    {
        stop("alignment.inputfile, umi.inputfile and gRNA.file
            are all required!")
    }
    if (length(alignment.inputfile) != length(umi.inputfile))
    {
       stop("number of alignment.inputfiles and number of
           umi.inputfiles should be the same")
    }
    if (length(alignment.inputfile) > 2 || length(umi.inputfile) > 2)
    {
       warning("Only the first 2 input files are used for
           the subsequent analysis")
    }
    n.files <- min(length(alignment.inputfile),
        length(umi.inputfile))
    n.files <- min(n.files, 2)
    if (class(gRNA.file) != "DNAStringSet")
        gRNAName <- gsub(".fa", "", basename(gRNA.file))
    else
        gRNAName <- names(gRNA.file)
    for (i in 1:length(gRNA.file))
    {
        if (is.na(gRNAName[i]) || nchar(gRNAName)[i]  == 0)
            gRNAName[i] <- paste("gRNAName", i, sep="")
    }
    cleavages.gr <- do.call(c, lapply(1:n.files, function(i)
    {
        cleavages <-
            getUniqueCleavageEvents(
            alignment.inputfile = alignment.inputfile[i],
            umi.inputfile = umi.inputfile[i],
            alignment.format = alignment.format,
            umi.header = umi.header, read.ID.col = read.ID.col,
            umi.col = umi.col, umi.sep = umi.sep,
            keep.chrM = keep.chrM,
            keep.R1only = keep.R1only, keep.R2only = keep.R2only,
            concordant.strand = concordant.strand,
            max.paired.distance = max.paired.distance,
            min.mapping.quality = min.mapping.quality,
            max.R1.len = max.R1.len, max.R2.len = max.R2.len,
            apply.both.max.len = apply.both.max.len,
            same.chromosome = same.chromosome,
            distance.inter.chrom = distance.inter.chrom,
            min.R1.mapped = min.R1.mapped,
            min.R2.mapped = min.R2.mapped,
            apply.both.min.mapped = apply.both.min.mapped,
            max.duplicate.distance = max.duplicate.distance,
            umi.plus.R1start.unique = umi.plus.R1start.unique,
            umi.plus.R2start.unique = umi.plus.R2start.unique,
            min.umi.count = min.umi.count,
            max.umi.count = max.umi.count,
            min.read.coverage = min.read.coverage,
            n.cores.max = n.cores.max)
        fileName <- gsub("bowtie2.", "",
            basename(alignment.inputfile[i]))
        fileName <- gsub("bowtie1.", "", fileName)
        fileName <- gsub("bowtie.", "", fileName)
        fileName <- gsub(".bed", "", fileName)
        fileName <- gsub(".sort", "", fileName)

        temp <- as.data.frame(cleavages$cleavage.gr)
        temp1 <- paste(temp[,1], temp[,5], temp[,2], sep = "")
        read.summary <- table(temp1)
        write.table(read.summary, file = paste(gRNAName, fileName,
            "ReadSummary.xls", sep = ""),
            sep = "\t", row.names = FALSE)
        seq.depth <- as.data.frame(table(cleavages$umi.count.summary$n))
        colnames(seq.depth)[1] <- c("UMIduplicateLevels")
        write.table(seq.depth, file = paste(gRNAName, fileName,
            "UMIsummary.xls", sep = ""),
            sep = "\t", row.names = FALSE)
	list(cleavages.gr = cleavages$cleavage.gr, read.summary = read.summary)
    }))
    message("Peak calling ...\n")
    if (n.files > 1)
         combined.gr <- c(cleavages.gr[[1]], cleavages.gr[[3]])
    else
         combined.gr <- cleavages.gr[[1]]
    peaks <- getPeaks(combined.gr, step = step,
        window.size = window.size, bg.window.size = bg.window.size,
        maxP = maxP, p.adjust.methods = p.adjust.methods,
        min.reads = min.reads, min.SNratio = min.SNratio)
    if (n.files >1)
    {
        peaks1 <- getPeaks(cleavages.gr[[1]], step = step,
            window.size = window.size, bg.window.size = bg.window.size,
            maxP = maxP, p.adjust.methods = p.adjust.methods,
            min.reads = min.reads.per.lib, min.SNratio = min.SNratio)
        peaks2 <- getPeaks(cleavages.gr[[3]], step = step,
            window.size = window.size, bg.window.size = bg.window.size,
            maxP = maxP, p.adjust.methods = p.adjust.methods,
            min.reads = min.reads.per.lib, min.SNratio = min.SNratio)
         #save(peaks1, file="peaks1.RData")
         #save(peaks2, file="peaks2.RData")

    }
    if (missing(outputDir))
    {
        outputDir <- getwd()
    }
    #write.table(as.data.frame(peaks$peaks),
    #    file = "testPeaks.xls", sep="\t", row.names=FALSE)
    #save(peaks, file="peaks.RData")
    message("combine plus and minus peaks ... \n")

    output.bedfile <- paste(gRNAName, "PlusMinusPeaksMerged.bed",
        sep = "-" )
    if (missing(outputDir) || outputDir == getwd())
    {
        outputDir <- paste(gRNAName,  "min", min.reads,
            "window", window.size, "step", step, "distance",
            distance.threshold, sep = "" )
    }
    if(!file.exists(outputDir))
        dir.create(outputDir)

    merged.gr<- mergePlusMinusPeaks(peaks.gr = peaks$peaks,
        distance.threshold = distance.threshold,
        max.overlap.plusSig.minusSig = max.overlap.plusSig.minusSig,
        output.bedfile = output.bedfile)
    append = FALSE
    if (length(merged.gr$mergedPeaks.gr) > 1)
    {
        write.table(cbind(name = names(merged.gr$mergedPeaks.gr),
            as.data.frame(merged.gr$mergedPeaks.gr)),
            file = file.path(outputDir, paste(gRNAName,
            "PlusMinusPeaksMerged.xls", sep = "-" )),
             sep="\t", quote = FALSE, row.names=FALSE, append = FALSE)
         append = TRUE
    }
   ##save(merged.gr, file="merged.gr.RData")
   #if (length(merged.gr) < 1)
       #merged.gr <- list(peaks.1strandOnly = peaks$peaks)
   message("keep peaks not in merged.gr but present in both peaks1 and peaks2\n")
   peaks.inboth1and2.gr <- GRanges()
    if (n.files >1 && length(peaks1$peaks) > 0 && length(peaks2$peaks) > 0)
    {
        cat("Find unmerged peaks with reads representation in both libraries\n")
        peaks.1strandOnly <- merged.gr$peaks.1strandOnly
        peaks.1strandOnly <-
            peaks.1strandOnly[peaks.1strandOnly$count >= (min.reads + 1)]
        if (length(names(peaks1$peaks)) < length(peaks1$peaks))
            names(peaks1$peaks) <- paste(seqnames(peaks1$peaks),
                start(peaks1$peaks), sep=":")
        if (length(names(peaks2$peaks)) < length(peaks2$peaks))
            names(peaks2$peaks) <- paste(seqnames(peaks2$peaks),
                start(peaks2$peaks), sep=":")
        peaks.in1.gr <- annotatePeakInBatch(peaks.1strandOnly, featureType = "TSS",
            AnnotationData = peaks1$peaks, output="overlap",
            PeakLocForDistance = "middle", FeatureLocForDistance = "middle",
            maxgap = 0)
        peaks.in1 <- unique(peaks.in1.gr[!is.na(peaks.in1.gr$feature)]$peak)
        peaks.in2.gr <- annotatePeakInBatch(peaks.1strandOnly, featureType = "TSS",
            AnnotationData = peaks2$peaks, output="overlap",
            PeakLocForDistance = "middle", FeatureLocForDistance = "middle",
            maxgap = 0)
        peaks.in2 <- unique(peaks.in2.gr[!is.na(peaks.in2.gr$feature)]$peak)

        peaks.inboth1and2 <- intersect(peaks.in1, peaks.in2)
        if (length(peaks.inboth1and2))
        {
            peaks.inboth1and2.gr <- peaks.1strandOnly[names(peaks.1strandOnly) %in%
                peaks.inboth1and2]
            bed.temp <- cbind(as.character(seqnames(peaks.inboth1and2.gr)),
                start(peaks.inboth1and2.gr),
                end(peaks.inboth1and2.gr), names(peaks.inboth1and2.gr),
                peaks.inboth1and2.gr$count, as.character(strand(peaks.inboth1and2.gr)))
            write.table(bed.temp, file = output.bedfile, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE, append = append)
            write.table(cbind(name = names(peaks.inboth1and2.gr),
                as.data.frame(peaks.inboth1and2.gr)),
                file = file.path(outputDir, paste(gRNAName,
                "PlusMinusPeaksMerged.xls", sep = "-" )),
                sep="\t", quote = FALSE, row.names=FALSE,
                col.names = FALSE, append = append)
            append = TRUE
         }
    }
    else if (!keepPeaksInBothStrandsOnly)
    {
       message("Find unmerged peaks with very high reads one-library protocol\n")
       peaks.1strandOnly.bed <- cbind(as.character(seqnames(merged.gr$peaks.1strandOnly)),
            start(merged.gr$peaks.1strandOnly),
            end(merged.gr$peaks.1strandOnly), names(merged.gr$peaks.1strandOnly),
            merged.gr$peaks.1strandOnly$count, as.character(strand(merged.gr$peaks.1strandOnly)))

       peaks.1strandOnly.bed <- peaks.1strandOnly.bed[as.numeric(as.character(
           peaks.1strandOnly.bed[,5])) >= min.peak.score.1strandOnly, ]
       write.table(peaks.1strandOnly.bed, file = output.bedfile, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE, append = append)
    }

    message("offtarget analysis ...\n")

    if(!file.exists(outputDir))
        dir.create(outputDir)

    offTargets <- offTargetAnalysisOfPeakRegions(gRNA = gRNA.file,
        peaks = output.bedfile,
        format = c(gRNA.format, "bed"),
        peaks.withHeader = FALSE, BSgenomeName = BSgenomeName,
        upstream = upstream, downstream = downstream,
        PAM.size = PAM.size, gRNA.size = gRNA.size,
        PAM =  PAM, PAM.pattern = PAM.pattern, max.mismatch = max.mismatch,
        outputDir = outputDir, overlap.gRNA.positions  = overlap.gRNA.positions,
        allowed.mismatch.PAM = allowed.mismatch.PAM, overwrite = overwrite,
        weights = weights,
        orderOfftargetsBy = keepTopOfftargetsBy,
        keepTopOfftargetsOnly = keepTopOfftargetsOnly,
        scoring.method = scoring.method,
        subPAM.activity = subPAM.activity,
        subPAM.position = subPAM.position,
        PAM.location = PAM.location,
        mismatch.activity.file = mismatch.activity.file
    )
    cat("Done with offtarget search!\n")
    offTargets <- subset(offTargets, !is.na(offTargets$offTarget))
    if (dim(offTargets)[1] == 0)
        stop("No offtargets found with the searching criteria!")
    cat("Add gene and exon information to offTargets ....\n")
    if (!missing(txdb) && (class(txdb) == "TxDb" ||
        class(txdb) == "TranscriptDb"))
    {
        offTargets <- annotateOffTargets(offTargets, txdb, orgAnn)
    }

    colnames(offTargets)[colnames(offTargets) == "n.mismatch"] <- "n.guide.mismatch"
    colnames(offTargets)[colnames(offTargets) == "name"] <- "gRNA.name"

    cat("Extract PAM sequence and n.PAM.mismatch. \n")
    if (PAM.location == "3prime")
    {
        PAM.sequence <- substr(offTargets$offTarget_sequence,
            gRNA.size + 1, gRNA.size + PAM.size)
    }
    else
    {
        PAM.sequence <- substr(offTargets$offTarget_sequence,
            1,  PAM.size)
    }

    n.PAM.mismatch <- unlist(lapply(DNAStringSet(PAM.sequence), function(i) {
        neditAt(i, DNAString(PAM), fixed=FALSE)
    }))
    exclude.col <- which(colnames(offTargets) %in%
        c("names", "targetSeqName", "peak_start", "peak_end", "peak_strand"))
    offTargets <- offTargets[, -exclude.col]
    ind1 <- which(colnames(offTargets) == "n.guide.mismatch")
    ind3 <- which(colnames(offTargets) == "mismatch.distance2PAM")
    ind.start1 <- min(ind1, ind3)
    ind2 <- ind.start1 + 1
    ind.start2 <- max(ind1, ind3)
    ind4 <- ind.start2 + 1
    if (ind.start2 > ind2)
    {
        offTargets <- cbind(offTargets[,1:ind.start1], n.PAM.mismatch = n.PAM.mismatch,
            offTargets[,ind2:ind.start2], PAM.sequence = PAM.sequence,
            offTargets[,ind4:dim(offTargets)[2]])
    }
    else
    {
        name.ind2 <- colnames(offTargets)[ind2]
        offTargets <- cbind(offTargets[,1:ind.start1], n.PAM.mismatch = n.PAM.mismatch,
            offTargets[,ind2], PAM.sequence = PAM.sequence,
            offTargets[,ind4:dim(offTargets)[2]])
        colnames(offTargets)[ind2+1] <- name.ind2
    }
    offTargets <- offTargets[order(offTargets[,which(colnames(offTargets) ==orderOfftargetsBy)],
        decreasing = descending), ]
    write.table(offTargets,
        file = file.path(outputDir,"offTargetsInPeakRegions.xls"),
        sep="\t", row.names = FALSE)

    message("Please check output file in directory ", outputDir , "\n")
    if (n.files > 1)
        list(offTargets = offTargets, merged.peaks = merged.gr$mergedPeaks.gr,
            peaks = peaks$peaks, uniqueCleavages = combined.gr,
            read.summary = list(s1 = cleavages.gr[[2]], s2 = cleavages.gr[[4]]),
            peaks.inboth1and2.gr = peaks.inboth1and2.gr)
    else
         list(offTargets = offTargets, merged.peaks = merged.gr$mergedPeaks.gr,
            peaks = peaks$peaks, uniqueCleavages = combined.gr,
            read.summary = cleavages.gr[[2]])
}

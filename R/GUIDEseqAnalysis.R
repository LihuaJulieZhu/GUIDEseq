GUIDEseqAnalysis <- function(alignment.inputfile,
    umi.inputfile,
    alignment.format = c("auto", "bam", "bed"),
    umi.header = FALSE,
    read.ID.col = 1, 
    umi.col = 2,
    umi.sep = "\t",
    BSgenomeName,
    gRNA.file,
    outputDir,
    n.cores.max = 6,
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
    window.size = 20L,
    step = 20L,
    bg.window.size = 5000L,
    min.reads = 5L,
    min.reads.per.lib = 1L,
    min.SNratio = 2,
    maxP = 0.05, 
    stats = c("poisson", "nbinom"),
    p.adjust.methods = 
        c( "none", "BH", "holm", "hochberg", "hommel", "bonferroni",
        "BY", "fdr"),
    distance.threshold = 40,
    plus.strand.start.gt.minus.strand.end = TRUE,
    gRNA.format = "fasta",
    overlap.gRNA.positions = c(17,18),
    upstream = 50,
    downstream = 50,
    PAM.size = 3,
    gRNA.size = 20,
    PAM = "NGG",
    PAM.pattern = "(NAG|NGG|NGA)$",
    max.mismatch = 6,
    allowed.mismatch.PAM = 2,
    overwrite = TRUE,
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,0.804, 0.685, 0.583),
    orderOfftargetsBy = c("predicted_cleavage_score", "n.mismatch"),
    descending = c(TRUE, FALSE),
    keepTopOfftargetsOnly = TRUE)
{
    alignment.format <- match.arg(alignment.format)
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
    gRNAName <- gsub(".fa", "", basename(gRNA.file))
    cleavages.gr <- do.call(c, lapply(1:n.files, function(i)
    {
        cleavages <-
            getUniqueCleavageEvents(
            alignment.inputfile = alignment.inputfile[i],
            umi.inputfile = umi.inputfile[i], 
            alignment.format = alignment.format,
            umi.header = umi.header, read.ID.col = read.ID.col,
            umi.col = umi.col, umi.sep = umi.sep,
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
    write.table(as.data.frame(peaks$peaks),
        file = file.path(outputDir, paste(gRNAName, 
        "peaks.xls", sep = "-" )), sep="\t", row.names=FALSE)

    message("combine plus and minus peaks ... \n")

    output.bedfile <- paste(gRNAName, "PlusMinusPeaksMerged.bed",
        sep = "-" )
    merged.gr<- mergePlusMinusPeaks(peaks.gr = peaks$peaks,
        distance.threshold = distance.threshold, step = step,
        output.bedfile = output.bedfile)
   #save(peaks, file="peaks.RData")
   #save(merged.gr, file="merged.gr.RData")
####### keep peaks not in merged.gr but present in both peaks1 and peaks2
    if (n.files >1)
    {
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
        peaks.inboth1and2.gr <- peaks.1strandOnly[names(peaks.1strandOnly) %in%
            peaks.inboth1and2]
        bed.temp <- cbind(as.character(seqnames(peaks.inboth1and2.gr)), 
            start(peaks.inboth1and2.gr),
            end(peaks.inboth1and2.gr), names(peaks.inboth1and2.gr),
            peaks.inboth1and2.gr$count, as.character(strand(peaks.inboth1and2.gr)))
       
        write.table(bed.temp, file = output.bedfile, sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
    }
    write.table(cbind(name = names(merged.gr$mergedPeaks.gr),
        as.data.frame(merged.gr$mergedPeaks.gr)),
        file = paste(gRNAName, "PlusMinusPeaksMerged.xls",
        sep = "-" ), sep="\t", quote = FALSE, row.names=FALSE)
    if (n.files > 1)
    {
        write.table(cbind(name = names(peaks.inboth1and2.gr),
            as.data.frame(peaks.inboth1and2.gr)),
            file = paste(gRNAName, "PlusMinusPeaksMerged.xls",
            sep = "-" ), sep="\t", quote = FALSE, row.names=FALSE, 
            col.names = FALSE, append = TRUE)
    } 
    message("offtarget analysis ...\n")
    inputFile1Path <- gRNA.file
    inputFile2Path <- output.bedfile

    if (missing(outputDir) || outputDir == getwd())
    {
        outputDir <- paste(gRNAName,  "min", min.reads,
            "window", window.size, "step", step, "distance",
            distance.threshold, sep = "" )
    }
    if(!file.exists(outputDir))
        dir.create(outputDir)
    write.table(cbind(name = names(merged.gr$mergedPeaks.gr),
        as.data.frame(merged.gr$mergedPeaks.gr)),
        file = file.path(outputDir, paste(gRNAName, 
        "PlusMinusPeaksMerged.xls", sep = "-" )),
        sep="\t", row.names=FALSE)

    offTargets <- offTargetAnalysisOfPeakRegions(gRNA = gRNA.file,
        peaks = output.bedfile,
        format = c(gRNA.format, "bed"),
        peaks.withHeader = FALSE, BSgenomeName = BSgenomeName,
        upstream = upstream, downstream = downstream,
        PAM.size = PAM.size, gRNA.size = gRNA.size,
        PAM =  PAM, PAM.pattern = PAM.pattern, max.mismatch = max.mismatch,
        outputDir = outputDir,overlap.gRNA.positions  = overlap.gRNA.positions,
        allowed.mismatch.PAM = allowed.mismatch.PAM, overwrite = overwrite,
        weights = weights,
        orderOfftargetsBy = orderOfftargetsBy,
        descending = descending,
        keepTopOfftargetsOnly = keepTopOfftargetsOnly
    )

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

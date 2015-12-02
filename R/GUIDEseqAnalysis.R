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
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,0.804, 0.685, 0.583))
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
    cleavages <-
        getUniqueCleavageEvents(alignment.inputfile = alignment.inputfile,
        umi.inputfile = umi.inputfile, alignment.format = alignment.format,
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
    gRNAName <- gsub(".fa", "", basename(gRNA.file))
    fileName <- gsub("bowtie2.", "", basename(alignment.inputfile))
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

    message("Peak calling ...\n")
    peaks <- getPeaks(cleavages$cleavage.gr, step = step,
        window.size = window.size, bg.window.size = bg.window.size,
       # n.cores.max = n.cores.max,
        maxP = maxP, p.adjust.methods = p.adjust.methods,
        min.reads = min.reads, min.SNratio = min.SNratio)

    if (missing(outputDir))
    {
        outputDir <- getwd()
    }
    write.table(as.data.frame(peaks$peaks),
        file = file.path(outputDir, paste(gRNAName, fileName,
        "peaks.xls", sep = "-" )), sep="\t", row.names=FALSE)

    message("combine plus and minus peaks ... \n")

    output.bedfile <- paste(gRNAName, fileName, "PlusMinusPeaksMerged.bed",
        sep = "-" )
    merged.gr<- mergePlusMinusPeaks(peaks.gr = peaks$peaks,
        distance.threshold = distance.threshold, step = step,
        output.bedfile = output.bedfile)

    write.table(cbind(name = names(merged.gr$mergedPeaks.gr),
        as.data.frame(merged.gr$mergedPeaks.gr)),
        file = paste(gRNAName, fileName, "PlusMinusPeaksMerged.xls",
        sep = "-" ), sep="\t", row.names=FALSE)

    message("offtarget analysis ...\n")
    inputFile1Path <- gRNA.file
    inputFile2Path <- output.bedfile

    if (missing(outputDir) || outputDir == getwd())
    {
        outputDir <- paste(gRNAName, fileName, "min", min.reads,
            "window", window.size, "step", step, "distance",
            distance.threshold, sep = "" )
    }
    if(!file.exists(outputDir))
        dir.create(outputDir)
    write.table(cbind(name = names(merged.gr$mergedPeaks.gr),
        as.data.frame(merged.gr$mergedPeaks.gr)),
        file = file.path(outputDir, paste(gRNAName, fileName,
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
        weights = weights
    )

    message("Please check output file in directory ", outputDir , "\n")
    list(offTargets = offTargets, merged.peaks = merged.gr$mergedPeaks.gr,
        peaks = peaks$peaks, uniqueCleavages = cleavages$cleavages.gr,
        read.summary = read.summary)
}

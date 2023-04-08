#' Analysis pipeline for GUIDE-seq dataset
#'
#' A wrapper function that uses the UMI sequence plus the first few bases of
#' each sequence from R1 reads to estimate the starting sequence library, piles
#' up reads with a user defined window and step size, identify the insertion
#' sites (proxy of cleavage sites), merge insertion sites from plus strand and
#' minus strand, followed by off target analysis of extended regions around the
#' identified insertion sites.
#'
#'
#' @param alignment.inputfile The alignment file. Currently supports bam and
#' bed output file with CIGAR information.  Suggest run the workflow
#' binReads.sh, which sequentially runs barcode binning, adaptor removal,
#' alignment to genome, alignment quality filtering, and bed file conversion.
#' Please download the workflow function and its helper scripts at
#' http://mccb.umassmed.edu/GUIDE-seq/binReads/
#' @param umi.inputfile A text file containing at least two columns, one is the
#' read identifier and the other is the UMI or UMI plus the first few bases of
#' R1 reads. Suggest use getUMI.sh to generate this file. Please download the
#' script and its helper scripts at http://mccb.umassmed.edu/GUIDE-seq/getUMI/
#' @param alignment.format The format of the alignment input file. Default bed
#' file format. Currently only bed file format is supported, which is generated
#' from binReads.sh
#' @param umi.header Indicates whether the umi input file contains a header
#' line or not. Default to FALSE
#' @param read.ID.col The index of the column containing the read identifier in
#' the umi input file, default to 1
#' @param umi.col The index of the column containing the umi or umi plus the
#' first few bases of sequence from the R1 reads, default to 2
#' @param umi.sep column separator in the umi input file, default to tab
#' @param BSgenomeName BSgenome object. Please refer to available.genomes in
#' BSgenome package. For example, BSgenome.Hsapiens.UCSC.hg19 for hg19,
#' BSgenome.Mmusculus.UCSC.mm10 for mm10, BSgenome.Celegans.UCSC.ce6 for ce6,
#' BSgenome.Rnorvegicus.UCSC.rn5 for rn5, BSgenome.Drerio.UCSC.danRer7 for Zv9,
#' and BSgenome.Dmelanogaster.UCSC.dm3 for dm3
#' @param gRNA.file gRNA input file path or a DNAStringSet object that contains
#' the target sequence (gRNA plus PAM)
#' @param outputDir the directory where the off target analysis and reports
#' will be written to
#' @param n.cores.max Indicating maximum number of cores to use in multi core
#' mode, i.e., parallel processing, default 1 to disable multicore processing
#' for small dataset.
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
#' @param min.umi.count To specify the minimum total count for a umi at the
#' genome level to be included in the subsequent analysis. For example, with
#' min.umi.count set to 2, if a umi only has 1 read in the entire genome, then
#' that umi will be excluded for the subsequent analysis.  Please adjust it to
#' a higher number for deeply sequenced library and vice versa.
#' @param max.umi.count To specify the maximum count for a umi to be included
#' in the subsequent analysis.  Please adjust it to a higher number for deeply
#' sequenced library and vice versa.
#' @param min.read.coverage To specify the minimum coverage for a read UMI
#' combination to be included in the subsequent analysis.  Please note that
#' this is different from min.umi.count which is less stringent.
#' @param apply.both.max.len Specify whether to apply maximum length
#' requirement to both R1 and R2 reads, default FALSE
#' @param same.chromosome Specify whether the paired reads are required to
#' align to the same chromosome, default TRUE
#' @param distance.inter.chrom Specify the distance value to assign to the
#' paired reads that are aligned to different chromosome, default -1
#' @param min.R1.mapped The minimum mapped R1 length to be considered for
#' downstream analysis, default 30 bp.
#' @param min.R2.mapped The minimum mapped R2 length to be considered for
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
#' @param window.size window size to calculate coverage
#' @param step step size to calculate coverage
#' @param bg.window.size window size to calculate local background
#' @param min.reads minimum number of reads to be considered as a peak
#' @param min.reads.per.lib minimum number of reads in each library (usually
#' two libraries) to be considered as a peak
#' @param min.peak.score.1strandOnly Specify the minimum number of reads for a
#' one-strand only peak to be included in the output. Applicable when set
#' keepPeaksInBothStrandsOnly to FALSE and there is only one library per sample
#' @param min.SNratio Specify the minimum signal noise ratio to be called as
#' peaks, which is the coverage normalized by local background.
#' @param maxP Specify the maximum p-value to be considered as significant
#' @param stats Statistical test, currently only poisson is implemented
#' @param p.adjust.methods Adjustment method for multiple comparisons, default
#' none
#' @param distance.threshold Specify the maximum gap allowed between the plus
#' strand and the negative strand peak, default 40. Suggest set it to twice of
#' window.size used for peak calling.
#' @param max.overlap.plusSig.minusSig Specify the cushion distance to allow
#' sequence error and inprecise integration Default to 30 to allow at most 10
#' (30-window.size 20) bp (half window) of minus-strand peaks on the right side
#' of plus-strand peaks. Only applicable if
#' plus.strand.start.gt.minus.strand.end is set to TRUE.
#' @param plus.strand.start.gt.minus.strand.end Specify whether plus strand
#' peak start greater than the paired negative strand peak end. Default to TRUE
#' @param keepPeaksInBothStrandsOnly Indicate whether only keep peaks present
#' in both strands as specified by plus.strand.start.gt.minus.strand.end,
#' max.overlap.plusSig.minusSig and distance.threshold.
#' @param gRNA.format Format of the gRNA input file. Currently, fasta is
#' supported
#' @param PAM.size PAM length, default 3
#' @param gRNA.size The size of the gRNA, default 20
#' @param PAM PAM sequence after the gRNA, default NGG
#' @param overlap.gRNA.positions The required overlap positions of gRNA and
#' restriction enzyme cut site, default 17 and 18 for SpCas9.
#' @param max.mismatch Maximum mismatch to the gRNA (not including mismatch to
#' the PAM) allowed in off target search, default 6
#' @param PAM.pattern Regular expression of protospacer-adjacent motif (PAM),
#' default NNN$. Alternatively set it to (NAG|NGG|NGA)$ for off target search
#' @param allowed.mismatch.PAM Maximum number of mismatches allowed for the PAM
#' sequence plus the number of degenerate sequence in the PAM sequence, default
#' to 2 for NGG PAM
#' @param upstream upstream offset from the peak start to search for off
#' targets, default 25 suggest set it to window size
#' @param downstream downstream offset from the peak end to search for off
#' targets, default 25 suggest set it to window size
#' @param overwrite overwrite the existing files in the output directory or
#' not, default FALSE
#' @param weights a numeric vector size of gRNA length, default c(0, 0, 0.014,
#' 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732,
#' 0.828, 0.615, 0.804, 0.685, 0.583) for SPcas9 system, which is used in Hsu
#' et al., 2013 cited in the reference section. Please make sure that the
#' number of elements in this vector is the same as the gRNA.size, e.g., pad 0s
#' at the beginning of the vector.
#' @param orderOfftargetsBy Criteria to order the offtargets, which works
#' together with the descending parameter
#' @param descending Indicate the output order of the offtargets, i.e., in the
#' descending or ascending order.
#' @param keepTopOfftargetsOnly Output all offtargets or the top offtarget
#' using the keepOfftargetsBy criteria, default to the top offtarget
#' @param keepTopOfftargetsBy Output the top offtarget for each called peak
#' using the keepTopOfftargetsBy criteria, If set to predicted_cleavage_score,
#' then the offtargets with the highest predicted cleavage score will be
#' retained If set to n.mismatch, then the offtarget with the lowest number of
#' mismatch to the target sequence will be retained
#' @param scoring.method Indicates which method to use for offtarget cleavage
#' rate estimation, currently two methods are supported, Hsu-Zhang and CFDscore
#' @param subPAM.activity Applicable only when scoring.method is set to
#' CFDscore A hash to represent the cleavage rate for each alternative sub PAM
#' sequence relative to preferred PAM sequence
#' @param subPAM.position Applicable only when scoring.method is set to
#' CFDscore The start and end positions of the sub PAM. Default to 22 and 23
#' for SP with 20bp gRNA and NGG as preferred PAM
#' @param PAM.location PAM location relative to gRNA. For example, default to
#' 3prime for spCas9 PAM.  Please set to 5prime for cpf1 PAM since it's PAM is
#' located on the 5 prime end
#' @param mismatch.activity.file Applicable only when scoring.method is set to
#' CFDscore A comma separated (csv) file containing the cleavage rates for all
#' possible types of single nucleotide mismatche at each position of the gRNA.
#' By default, use the supplemental Table 19 from Doench et al., Nature
#' Biotechnology 2016
#' @param bulge.activity.file Used for predicting indel effect on offtarget
#' cleavage score. An excel file with the second sheet for deletion activity
#' and the third sheet for Insertion. By default, use the supplemental Table 19
#' from Doench et al., Nature Biotechnology 2016
#' @param txdb TxDb object, for creating and using TxDb object, please refer to
#' GenomicFeatures package. For a list of existing TxDb object, please search
#' for annotation package starting with Txdb at
#' http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData,
#' such as TxDb.Rnorvegicus.UCSC.rn5.refGene for rat,
#' TxDb.Mmusculus.UCSC.mm10.knownGene for mouse,
#' TxDb.Hsapiens.UCSC.hg19.knownGene for human,
#' TxDb.Dmelanogaster.UCSC.dm3.ensGene for Drosophila and
#' TxDb.Celegans.UCSC.ce6.ensGene for C.elegans
#' @param orgAnn organism annotation mapping such as org.Hs.egSYMBOL in
#' org.Hs.eg.db package for human
#' @param mat nucleotide substitution matrix. Function
#' nucleotideSubstitutionMatrix can be used for creating customized nucleotide
#' substitution matrix. By default, match = 1, mismatch = -1, and baseOnly =
#' TRUE Only applicalbe with includeBulge set to TRUE
#' @param includeBulge indicates whether including offtargets with indels
#' default to FALSE
#' @param max.n.bulge offtargets with at most this number of indels to be
#' included in the offtarget list. Only applicalbe with includeBulge set to
#' TRUE
#' @param min.peak.score.bulge default to 60. Set it to a higher number to
#' speed up the alignment with bulges. Any peaks with 
#' peak.score less than min.peak.score.bulge will not be included in
#' the offtarget analysis with bulges. However, all peaks
#' will be included in the offtarget analysis with mismatches. 
#' @param removeDuplicate default to TRUE. Set it to FALSE if PCR duplicates
#' not to be removed for testing purpose
#' @param resume default to FALSE to restart the analysis. set it TRUE to
#' resume an analysis.
#' @param ignoreTagmSite default to FALSE. To collapse
#' reads with the same integration site and UMI but with different 
#' tagmentation site, set the option to TRUE.
#' @param ignoreUMI default to FALSE. To collapse reads with the same 
#' integration and tagmentation site but with different UMIs, 
#' set the option to TRUE and retain the UMI that appears most frequently
#'  for each combination of integration and tagmentation site. 
#'  In case of ties, randomly select one UMI.
#'
#' @return \item{offTargets}{ a data frame, containing all input peaks with
#' potential gRNA binding sites, mismatch number and positions, alignment to
#' the input gRNA and predicted cleavage score.}
#' \item{merged.peaks}{merged
#' peaks as GRanges with count (peak height), bg (local background), SNratio
#' (signal noise ratio), p-value, and option adjusted p-value }
#' \item{peaks}{GRanges with count (peak height), bg (local background), SNratio (signal
#' noise ratio), p-value, and option adjusted p-value }
#' \item{uniqueCleavages}{Cleavage sites with one site per UMI as GRanges with
#' metadata column total set to 1 for each range}
#' \item{read.summary}{One table
#' per input mapping file that contains the number of reads for each chromosome
#' location}
#' \item{sequence.depth}{sequence depth in the input alignment files}
#' @author Lihua Julie Zhu
#' @seealso getPeaks
#' @references Lihua Julie Zhu, Michael Lawrence, Ankit Gupta, Herve Pages, Alper Ku- cukural, Manuel Garber and Scot A. Wolfe. GUIDEseq: a bioconductor package to analyze GUIDE-Seq datasets for CRISPR-Cas nucleases. BMC Genomics. 2017. 18:379
#' @keywords misc
#' @examples
#'
#' if(interactive())
#'     {
#'         library("BSgenome.Hsapiens.UCSC.hg19")
#'         umiFile <- system.file("extdata", "UMI-HEK293_site4_chr13.txt",
#'            package = "GUIDEseq")
#'         alignFile <- system.file("extdata","bowtie2.HEK293_site4_chr13.sort.bam" ,
#'             package = "GUIDEseq")
#'         gRNA.file <- system.file("extdata","gRNA.fa", package = "GUIDEseq")
#'         guideSeqRes <- GUIDEseqAnalysis(
#'             alignment.inputfile = alignFile,
#'             umi.inputfile = umiFile, gRNA.file = gRNA.file,
#'             orderOfftargetsBy = "peak_score",
#'             descending = TRUE,
#'             keepTopOfftargetsBy = "predicted_cleavage_score",
#'             scoring.method = "CFDscore",
#'             BSgenomeName = Hsapiens, min.reads = 80, n.cores.max = 1)
#'         guideSeqRes$offTargets
#'         names(guideSeqRes)
#'    }
#'
#' @importFrom tidyr unnest
#' @importFrom utils read.table write.table
#' @importFrom dplyr mutate group_by add_count select
#' filter '%>%' rename ungroup n
#' @importFrom GenomicRanges GRanges start end
#' strand seqnames
#' @importFrom Biostrings DNAString DNAStringSet neditAt
#' readDNAStringSet
#' @importFrom ChIPpeakAnno annotatePeakInBatch
#' @importFrom hash hash
#'
#' @export GUIDEseqAnalysis

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
    min.umi.count = 1L,
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
    upstream = 25L,
    downstream = 25L,
    PAM.size = 3L,
    gRNA.size = 20L,
    PAM = "NGG",
    PAM.pattern = "NNN$",
    max.mismatch = 6L,
    allowed.mismatch.PAM = 2L,
    overwrite = TRUE,
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,0.804, 0.685, 0.583),
    orderOfftargetsBy = c("peak_score", "predicted_cleavage_score", "n.guide.mismatch"),
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
     bulge.activity.file = system.file("extdata",
         "NatureBiot2016SuppTable19DoenchRoot.xlsx",
                  package = "GUIDEseq"),
     txdb, orgAnn,
     mat,
     includeBulge = FALSE,
     max.n.bulge = 2L,
     min.peak.score.bulge = 60L,
     removeDuplicate = TRUE,
     resume = FALSE,
     ignoreTagmSite = FALSE,
     ignoreUMI = FALSE)
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
    {
        #gRNAName <- gsub(".fa", "", basename(gRNA.file))
        gRNA <- readDNAStringSet(gRNA.file)
    }
    else
    {
        gRNA <- gRNA.file
    }

    gRNAName <- names(gRNA)
    if (PAM.location == "3prime")
        gRNA <- as.character(substr(as.character(gRNA), 1, gRNA.size))
    else if (nchar(as.character(gRNA)) == gRNA.size + PAM.size)
        gRNA <- as.character(substr(as.character(gRNA), PAM.size + 1,
                       gRNA.size + PAM.size))
    else if (nchar(as.character(gRNA)) != gRNA.size)
        stop("gRNA should include protospace with/without PAM")

    for (i in 1:length(gRNA.file))
    {
        if (length(gRNAName[i]) == 0)
            gRNAName[i] <- paste("gRNAName", i, sep="")
    }
    if (missing(outputDir) || outputDir == getwd())
    {
        outputDir <- paste(gsub(".sorted", "", gsub(".bam", "",
                                 basename(alignment.inputfile[1]))),
                            paste(gRNAName,  "min", min.reads,
                           "window", window.size, "step", step, "distance",
                           distance.threshold, sep = "" ), sep ="_")
    }
    if(!file.exists(outputDir))
        dir.create(outputDir, recursive=TRUE)

    sampleName <- gsub(".", "",
                       gsub("minus", "", gsub("plus", "", gsub(".bam", "",
                       gsub("sort", "",
                       gsub("sorted", "",
                            gsub("bowtie", "",
                                gsub("bowtie1", "",
                                    gsub("bowtie2", "",
                        basename(alignment.inputfile[1]))))))))),
                      fixed = TRUE)

    output.bedfile <- file.path(outputDir, paste(sampleName,
                            paste(gRNAName, "PlusMinusPeaksMerged.bed",
                                  sep = "_" ), sep ="_"))
    skipPeakCalling <- TRUE
    if (!resume || !file.exists(output.bedfile)) {
         skipPeakCalling <- FALSE
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
               n.cores.max = n.cores.max,
               outputDir = outputDir,
               removeDuplicate = removeDuplicate,
               ignoreTagmSite = ignoreTagmSite,
               ignoreUMI = ignoreUMI)
            fileName <- gsub("bowtie2.", "",
               basename(alignment.inputfile[i]))
            fileName <- gsub("bowtie1.", "", fileName)
            fileName <- gsub("bowtie.", "", fileName)
            fileName <- gsub(".bed", "", fileName)
            fileName <- gsub(".sorted", "", fileName)
            fileName <- gsub(".bam", "", fileName)
            fileName <- gsub(".sort", "", fileName)

            temp <- as.data.frame(cleavages$cleavage.gr)
            temp1 <- paste(temp[,1], temp[,5], temp[,2], sep = "")
            read.summary <- table(temp1)
            write.table(read.summary,
                    file = file.path(outputDir,
                                     paste(gRNAName, fileName,
                    "ReadSummary.xls", sep = "")),
                  sep = "\t", row.names = FALSE)
            seq.depth <- as.data.frame(table(cleavages$umi.count.summary$n))
            colnames(seq.depth)[1] <- c("UMIduplicateLevels")
            write.table(seq.depth,
                    file = file.path(outputDir,
                                     paste(gRNAName, fileName,
                    "UMIsummary.xls", sep = "")),
                    sep = "\t", row.names = FALSE)
	     list(cleavages.gr = cleavages$cleavage.gr,
	         read.summary = read.summary,
	         sequence.depth = cleavages$sequence.depth)
         }))
        saveRDS(cleavages.gr, file = file.path(outputDir,"cleavages.RDS"))
        message("Peak calling ...\n")
        if (n.files > 1)
              combined.gr <- c(cleavages.gr[[1]], cleavages.gr[[4]])
        else
              combined.gr <- cleavages.gr[[1]]
        saveRDS(combined.gr, file.path(outputDir,"combined.gr.RDS"))
        
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
            peaks2 <- getPeaks(cleavages.gr[[4]], step = step,
                window.size = window.size, bg.window.size = bg.window.size,
                maxP = maxP, p.adjust.methods = p.adjust.methods,
                min.reads = min.reads.per.lib, min.SNratio = min.SNratio)
         #save(peaks1, file="peaks1.RData")
         #save(peaks2, file="peaks2.RData")
        }

    #write.table(as.data.frame(peaks$peaks),
    #    file = "testPeaks.xls", sep="\t", row.names=FALSE)
    #save(peaks, file="peaks.RData")
       message("combine plus and minus peaks ... \n")

       merged.gr<- mergePlusMinusPeaks(peaks.gr = peaks$peaks,
           distance.threshold = distance.threshold,
           max.overlap.plusSig.minusSig = max.overlap.plusSig.minusSig,
           output.bedfile = output.bedfile)
       append = FALSE
       if (length(merged.gr$mergedPeaks.gr) >= 1 && 
           class(merged.gr$mergedPeaks.gr) == "GRanges")
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
                write.table(bed.temp,
                        file = output.bedfile,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE,
                        quote = FALSE, append = append)
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
            write.table(peaks.1strandOnly.bed,
                   file = output.bedfile,
                   sep = "\t",
                   row.names = FALSE,
                   col.names = FALSE,
                   quote = FALSE, append = append)
         }
    } # if merged bed file does not exist
    message("offtarget analysis ...\n")
    offTargets <-  read.table(
            file = output.bedfile,
            sep = "\t",
            header = FALSE)
    skipOTAWithNoBulge <- TRUE
    if (!resume || !file.exists(file.path(outputDir,"offTargetsWithNoBulge.RDS")))
    {
       skipOTAWithNoBulge <- FALSE
       tryCatch(offTargets <- offTargetAnalysisOfPeakRegions(gRNA = gRNA.file,
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
            mismatch.activity.file = mismatch.activity.file),
        error = function(e) {
             message(e)
             return(NA)
        })

       if (length(which(colnames(offTargets) == "n.mismatch")) > 0)
       {
            colnames(offTargets)[colnames(offTargets) == "n.mismatch"] <-
                "n.guide.mismatch"
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

            n.PAM.mismatch <- unlist(lapply(DNAStringSet(PAM.sequence),
                                        function(i) {
                neditAt(i, DNAString(PAM), fixed=FALSE)
            }))
            exclude.col <- which(colnames(offTargets) %in%
                                 c("names", "targetSeqName",
                                   "peak_start", "peak_end", "peak_strand"))
            offTargets <- offTargets[, -exclude.col]
            ind1 <- which(colnames(offTargets) == "n.guide.mismatch")
            ind3 <- which(colnames(offTargets) == "mismatch.distance2PAM")
            ind.start1 <- min(ind1, ind3)
            ind2 <- ind.start1 + 1
            ind.start2 <- max(ind1, ind3)
            ind4 <- ind.start2 + 1
            if (ind.start2 > ind2)
            {
                offTargets <- cbind(offTargets[,1:ind.start1],
                                n.PAM.mismatch = n.PAM.mismatch,
                                offTargets[,ind2:ind.start2],
                                PAM.sequence = PAM.sequence,
                                offTargets[,ind4:dim(offTargets)[2]])
            }
            else
            {
                name.ind2 <- colnames(offTargets)[ind2]
                offTargets <- cbind(offTargets[,1:ind.start1],
                                n.PAM.mismatch = n.PAM.mismatch,
                                offTargets[,ind2],
                                PAM.sequence = PAM.sequence,
                                offTargets[,ind4:dim(offTargets)[2]])
                colnames(offTargets)[ind2+1] <- name.ind2
            }
        }
        saveRDS(offTargets, file = file.path(outputDir,"offTargetsWithNoBulge.RDS"))
    } # if offtarget analysis has not been performed previously
    else
    {
        offTargets <- readRDS(file.path(outputDir,"offTargetsWithNoBulge.RDS"))
        cleavages.gr <- readRDS(file.path(outputDir,"cleavages.RDS"))
        combined.gr <- readRDS(file.path(outputDir,"combined.gr.RDS"))
    }
    if (includeBulge)
    {
         if (length(which(colnames(offTargets) == "n.guide.mismatch")) > 0)
        {
            offTargets <- offTargets %>%
                select(offTarget,
                       peak_score,
                       predicted_cleavage_score,
                       gRNA.name,
                       gRNAPlusPAM,
                       offTarget_sequence,
                       guideAlignment2OffTarget,
                       offTargetStrand,
                       mismatch.distance2PAM,
                       n.PAM.mismatch,
                       n.guide.mismatch,
                       PAM.sequence,
                       offTarget_Start,
                       offTarget_End,
                       chromosome) %>%
                mutate(RNA.bulge = "",
                       DNA.bulge = "",
                       pos.RNA.bulge = "",
                       pos.DNA.bulge = "",
                       n.RNA.bulge = 0,
                       n.DNA.bulge = 0,
                       total.mismatch.bulge =
                           n.guide.mismatch + n.PAM.mismatch)
        }
        cat("Finding offtargets with bulges ...")
        bed.for.bulge <- read.table(output.bedfile, sep = "\t", 
             header = FALSE)
        bed.for.bulge <- bed.for.bulge[bed.for.bulge[,5] >= min.peak.score.bulge,]
        bulge.bedfile <- file.path(outputDir, paste(sampleName,
                            paste(gRNAName, "peaksForBulgeSearch.bed",
                                  sep = "_" ), sep ="_"))
        write.table(bed.for.bulge, file = bulge.bedfile,
             sep = "\t", col.names = FALSE, row.names = FALSE) 
        offTargets.bulge <- do.call(rbind, lapply(1:length(gRNAName), function(i) {
            temp <- offTargetAnalysisWithBulge(gRNA = gRNA[i], gRNA.name = gRNAName[i],
                                               peaks = bulge.bedfile,
                                               BSgenomeName = BSgenomeName,
                                               mismatch.activity.file = bulge.activity.file,
                                               peaks.withHeader = FALSE,
                                               PAM.size = PAM.size, gRNA.size = gRNA.size,
                                               PAM =  PAM, PAM.pattern = PAM.pattern,
                                               PAM.location = PAM.location,
                                               max.mismatch = max.mismatch,
                                               allowed.mismatch.PAM = allowed.mismatch.PAM,
                                               max.DNA.bulge = max.n.bulge,
                                               upstream = upstream,
                                               downstream = downstream)
            temp$score.bulge
        }))
        if (class(offTargets.bulge) == "data.frame" &&
            nrow(offTargets.bulge) > 0)
        {
           colnames(offTargets.bulge)[colnames(offTargets.bulge) ==
                                 "n.mismatch"] <- "n.guide.mismatch"
            offTargets.bulge <- offTargets.bulge[
                unlist(offTargets.bulge$n.deletion) > 0 |
                unlist(offTargets.bulge$n.insertion) > 0,]
            if (nrow(offTargets.bulge) > 0)
            {
                offTargets.b <- offTargets.bulge %>%
                  select(offTarget,
                   peak_score,
                   predicted_cleavage_score,
                   gRNA.name,
                   gRNAPlusPAM,
                   offTarget_sequence,
                   guideAlignment2OffTarget,
                   offTargetStrand,
                   mismatch.distance2PAM,
                   n.PAM.mismatch,
                   n.guide.mismatch,
                   PAM.sequence,
                   offTarget_Start,
                   offTarget_End,
                   chromosome,
                   gRNA.insertion, gRNA.deletion,
                   pos.insertion, pos.deletion,
                   n.insertion, n.deletion) %>%
                rename(RNA.bulge = gRNA.insertion,
                       DNA.bulge = gRNA.deletion,
                       pos.RNA.bulge  = pos.insertion,
                       pos.DNA.bulge = pos.deletion,
                       n.RNA.bulge = n.insertion,
                       n.DNA.bulge = n.deletion) %>%
               mutate(total.mismatch.bulge =
                       as.numeric(n.guide.mismatch) +
                          as.numeric(n.PAM.mismatch) +
                          as.numeric(n.RNA.bulge) +
                          as.numeric(n.DNA.bulge))

            offTargets.b.bk  <- do.call(cbind, lapply(1:ncol(offTargets.b),
                                                      function(i) {
                temp <- offTargets.b[,i]
                if (mode(temp[1]) == "list")
                {
                    unlist(lapply(1:length(temp), function(j) {
                     temp2 <- unlist(temp[j])
                        ifelse(length(temp2) > 1,
                              paste(temp2, collapse = ","), temp2)
                    }))
                }
                else
                {
                    temp
                }
            }))
            colnames(offTargets.b.bk) <- colnames(offTargets.b)
            offTargets.b <- offTargets.b.bk
            rm(offTargets.b.bk)

            offTargets.b <- as.data.frame(offTargets.b)
            offTargets.b[,which(colnames(offTargets.b) == "offTarget_End")] <-
                as.numeric(offTargets.b[,which(colnames(offTargets.b) ==
                                               "offTarget_End")] )
            offTargets.b[,which(colnames(offTargets.b) == "offTarget_Start")] <-
                as.numeric(offTargets.b[,which(colnames(offTargets.b) ==
                                               "offTarget_Start")] )

            offTargets.b$pos.RNA.bulge[is.na(offTargets.b$pos.RNA.bulge)] <- ""
            offTargets.b$pos.DNA.bulge[is.na(offTargets.b$pos.DNA.bulge)] <- ""
            # convert U back to T for visualization (U is for CFD score calculation)
            offTargets.b$DNA.bulge <- gsub("U", "T", offTargets.b$DNA.bulge)

            if (length(which(colnames(offTargets) == "n.guide.mismatch")) > 0)
            {
                offTargets <- rbind(offTargets.b, offTargets)
                offTargets <- subset(offTargets, 
                                    (as.numeric(offTargets$n.DNA.bulge) + 
                                        as.numeric(offTargets$n.RNA.bulge)) <= max.n.bulge)
            }
            else
            {
                offTargets <- offTargets.b
            }
          }
        }
    }
    cat("Done with offtarget search!\n")

    offTargets <- subset(offTargets, !is.na(offTargets$offTarget))
    
    if (nrow(offTargets) == 0)
    {
        message("No offtargets found with the searching criteria!")
    }
    else
    {
        cat("Add gene and exon information to offTargets ....\n")
        if (!missing(txdb) && (class(txdb) == "TxDb" ||
            class(txdb) == "TranscriptDb"))
        {
            offTargets <- annotateOffTargets(offTargets, txdb, orgAnn)
        }

        cat("Order offtargets. \n")

        orderOfftargetsBy <- intersect(orderOfftargetsBy,
                                       colnames(offTargets))
        if(length(orderOfftargetsBy) > 0)
            offTargets <- offTargets[order(offTargets[,which(
                colnames(offTargets) == orderOfftargetsBy)],
                    decreasing = descending), ]

        message("Add sequence depth information. \n")
        if (n.files == 1)
            offTargets <- offTargets %>%
                mutate(sequence.depth = cleavages.gr[[3]])
        else if (n.files == 2)
            offTargets <- offTargets %>%
                mutate(sequence.depth =
                           cleavages.gr[[3]] + cleavages.gr[[6]])
        message("Get the number of unique UMIs for each offtarget. \n")
        cnames <- colnames(offTargets)
        x <- makeGRangesFromDataFrame(offTargets,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=FALSE,
                                      seqinfo=NULL,
                                      seqnames.field="chromosome",
                                      start.field = "offTarget_Start",
                                      end.field = "offTarget_End",
                                      strand.field = "offTargetStrand",
                                      starts.in.df.are.0based = FALSE)

        names(combined.gr) <- paste0(combined.gr$umi,
                                     seq_along(1:length(combined.gr$umi)))
        ann.offtargets.gr <- annotatePeakInBatch(x, AnnotationData = combined.gr,
                                                 output = "overlap",
                                                 maxgap = max(upstream, downstream))
        offTargets <- as.data.frame(ann.offtargets.gr) %>%
            group_by(offTarget) %>%
            mutate(n.distinct.UMIs = length(unique(substr(feature,
                                           1, nchar(combined.gr$umi[1])))),
                   peak_score = n()) %>%
            rename(offTarget_Start = start,
                   offTarget_End = end,
                   offTargetStrand = strand,
                  chromosome = seqnames) %>%
                  ungroup  %>%
            select(!!cnames, n.distinct.UMIs) %>%
            unique %>% as.data.frame
        cat("Save offtargets. \n")
        write.table(offTargets, file =
            file.path(outputDir,"offTargetsInPeakRegions.xls"),
            sep = "\t", row.names = FALSE)
    }

    message("Please check output file in directory ", outputDir , "\n")

    if (skipPeakCalling)
    {
          list(offTargets = offTargets)
    }
    else if (n.files > 1)
    {
        list(offTargets = offTargets, merged.peaks = merged.gr$mergedPeaks.gr,
            peaks = peaks$peaks, uniqueCleavages = combined.gr,
            read.summary = list(s1 = cleavages.gr[[2]], s2 = cleavages.gr[[5]]),
            peaks.inboth1and2.gr = peaks.inboth1and2.gr,
            sequence.depth = list(s1 = cleavages.gr[[3]],
                                  s2 = cleavages.gr[[6]])
        )
    }
    else
    {
         list(offTargets = offTargets, merged.peaks = merged.gr$mergedPeaks.gr,
            peaks = peaks$peaks, uniqueCleavages = combined.gr,
            read.summary = cleavages.gr[[2]],
            sequence.depth = cleavages.gr[[3]])
    }
}

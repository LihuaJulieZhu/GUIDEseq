#' Analysis pipeline for PEtag-seq dataset
#'
#' A wrapper function that uses the UMI sequence plus the first few bases of
#' each sequence from R1 reads to estimate the starting sequence library, piles
#' up reads with a user defined window and step size, identify the insertion
#' sites (proxy of cleavage sites), merge insertion sites from plus strand and
#' minus strand, followed by off target analysis of extended regions around the
#' identified insertion sites. Detailed information on additional parameters
#' can be found in GUIDEseqAnalysis manual with help(GUIDEseqAnalysis).
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
#' @param BSgenomeName BSgenome object. Please refer to available.genomes in
#' BSgenome package. For example, BSgenome.Hsapiens.UCSC.hg19 for hg19,
#' BSgenome.Mmusculus.UCSC.mm10 for mm10, BSgenome.Celegans.UCSC.ce6 for ce6,
#' BSgenome.Rnorvegicus.UCSC.rn5 for rn5, BSgenome.Drerio.UCSC.danRer7 for Zv9,
#' and BSgenome.Dmelanogaster.UCSC.dm3 for dm3
#' @param gRNA.file gRNA input file path or a DNAStringSet object that contains
#' the target sequence (gRNA plus PAM)
#' @param outputDir the directory where the off target analysis and reports
#' will be written to
#' @param keepPeaksInBothStrandsOnly Indicate whether only keep peaks present
#' in both strands as specified by plus.strand.start.gt.minus.strand.end,
#' max.overlap.plusSig.minusSig and distance.threshold. Please see
#' GUIDEseqAnalysis for details of additional parameters. Default to FALSE for
#' any in vitro system, which needs to be set to TRUE for any in vivo system.
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
#' @param PAM.size PAM length, default 3
#' @param gRNA.size The size of the gRNA, default 20
#' @param overlap.gRNA.positions The required overlap positions of gRNA and
#' restriction enzyme cut site, default 17 and 18 for SpCas9.
#' @param PAM.location PAM location relative to gRNA. For example, default to
#' 3prime for spCas9 PAM.  Please set to 5prime for cpf1 PAM since it's PAM is
#' located on the 5 prime end
#' @param PBS.len Primer binding sequence length, default to 10.
#' @param HA.len Homology arm sequence length, default to 7.
#' @param ...  Any parameters in GUIDEseqAnalysis can be used for this
#' function. Please type help(GUIDEseqAnalysis for detailed information.
#' @return \item{offTargets}{ a data frame, containing all input peaks with
#' potential gRNA binding sites, mismatch number and positions, alignment to
#' the input gRNA, predicted cleavage score, PBS (primer binding sequence), and
#' HAseq (homology arm sequence).} \item{merged.peaks}{merged peaks as GRanges
#' with count (peak height), bg (local background), SNratio (signal noise
#' ratio), p-value, and option adjusted p-value } \item{peaks }{GRanges with
#' count (peak height), bg (local background), SNratio (signal noise ratio),
#' p-value, and option adjusted p-value } \item{uniqueCleavages}{Cleavage sites
#' with one site per UMI as GRanges with metadata column total set to 1 for
#' each range} \item{read.summary}{One table per input mapping file that
#' contains the number of reads for each chromosome location}
#' @author Lihua Julie Zhu
#' @seealso GUIDEseqAnalysis
#' @references Lihua Julie Zhu, Michael Lawrence, Ankit Gupta, Herve Pages, Alper Ku- cukural, Manuel Garber and Scot A. Wolfe. GUIDEseq: a bioconductor package to analyze GUIDE-Seq datasets for CRISPR-Cas nucleases. BMC Genomics. 2017. 18:379
#' @keywords misc
#' @examples
#'
#' if(!interactive())
#'     {
#'         library("BSgenome.Hsapiens.UCSC.hg19")
#'         library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'         library(org.Hs.eg.db)
#'         umiFile <- system.file("extdata", "UMI-HEK293_site4_chr13.txt",
#'            package = "GUIDEseq")
#'         alignFile <- system.file("extdata","bowtie2.HEK293_site4_chr13.sort.bam" ,
#'             package = "GUIDEseq")
#'         gRNA.file <- system.file("extdata","gRNA.fa", package = "GUIDEseq")
#'         PET.res <- PEtagAnalysis(
#'             alignment.inputfile = alignFile,
#'             umi.inputfile = umiFile,
#'             gRNA.file = gRNA.file,
#'             orderOfftargetsBy = "peak_score",
#'             descending = TRUE,
#'             keepTopOfftargetsBy = "predicted_cleavage_score",
#'             scoring.method = "CFDscore",
#'             BSgenomeName = Hsapiens,
#'             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'             orgAnn = org.Hs.egSYMBOL,
#'             outputDir = "PEtagTestResults",
#'             min.reads = 80, n.cores.max = 1,
#'             keepPeaksInBothStrandsOnly = FALSE,
#'             PBS.len = 10L,
#'             HA.len = 7L
#'             )
#'         PET.res$offTargets
#'         names(PET.res)
#'    }
#'
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom utils write.table
#' @importFrom BSgenome getSeq
#'
#' @export PEtagAnalysis
PEtagAnalysis <- function(alignment.inputfile,
    umi.inputfile,
    BSgenomeName,
    gRNA.file,
    outputDir,
    keepPeaksInBothStrandsOnly = FALSE,
    txdb,
    orgAnn,
    PAM.size = 3L,
    gRNA.size = 20L,
    overlap.gRNA.positions = c(17,18),
    PAM.location = "3prime",
    PBS.len = 10L,
    HA.len = 7L,
    ...)
{
   res <- GUIDEseqAnalysis(alignment.inputfile = alignment.inputfile,
               umi.inputfile = umi.inputfile,
               BSgenomeName = BSgenomeName,
               gRNA.file = gRNA.file,
               keepPeaksInBothStrandsOnly = keepPeaksInBothStrandsOnly,
               txdb = txdb,
               orgAnn = orgAnn,
               PAM.size = PAM.size,
               gRNA.size = gRNA.size,
               overlap.gRNA.positions = overlap.gRNA.positions,
               outputDir = outputDir,
               PAM.location = PAM.location,
               ...)
   if (PAM.location == "3prime")
   {
      offtargets <- res[[1]]
      baseBeforegRNA <- PBS.len - overlap.gRNA.positions[1]
      baseAfterPAM <- HA.len - PAM.size - gRNA.size + overlap.gRNA.positions[1]

      if (dim(offtargets)[1] > 0)
      {
        chr <- as.character(offtargets$chromosome)
        strand <- as.character(offtargets$offTargetStrand)
        Start <- ifelse(strand=="-",
              as.numeric(as.character( offtargets$offTarget_Start)) - baseAfterPAM,
              as.numeric(as.character( offtargets$offTarget_Start)) - baseBeforegRNA)
        End <- ifelse(strand=="-",
              as.numeric(as.character( offtargets$offTarget_End)) + as.numeric(baseBeforegRNA),
              as.numeric(as.character( offtargets$offTarget_End)) + as.numeric(baseAfterPAM))
       }
       starts <- unlist(apply(cbind(Start,1), 1, max))
       ends <- unlist(apply(cbind(End, seqlengths(BSgenomeName)[chr]), 1,min))
       extendedSequence <- getSeq(BSgenomeName, names = chr, start = starts,
             end = ends, strand = strand, width = NA, as.character = TRUE)
       PBS <- substr(extendedSequence, 1,  PBS.len)
       HA <- substr(extendedSequence, PBS.len + 1, PBS.len + HA.len)
       offtargets <- cbind(offtargets, PBS = PBS, HAseq = HA)
       res[[1]] <- offtargets
       write.table(offtargets,
          file = file.path(outputDir,"offTargetsInPeakRegions.xls"),
          sep="\t", row.names = FALSE)
  }
  res
}



#' Analysis of GUIDE-seq
#' 
#' The package includes functions to retain one read per unique molecular
#' identifier (UMI), filter reads lacking integration oligo sequence, identify
#' peak locations (cleavage sites) and heights, merge peaks, perform off-target
#' search using the input gRNA. This package leverages CRISPRseek and
#' ChIPpeakAnno packages.
#' 
#' \tabular{ll}{ Package: \tab GUIDEseq\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2015-09-04\cr License: \tab GPL (>= 2)\cr } Function
#' GUIDEseqAnalysis integrates all steps of GUIDE-seq analysis into one
#' function call
#' 
#' @name GUIDEseq-package
#' @aliases GUIDEseq-package GUIDEseq
#' @docType package
#' @author Lihua Julie Zhu Maintainer:julie.zhu@@umassmed.edu
#' @seealso GUIDEseqAnalysis
#' @references Shengdar Q Tsai and J Keith Joung et al. GUIDE-seq enables
#' genome-wide profiling of off-target cleavage by CRISPR-Cas nucleases. Nature
#' Biotechnology 33, 187 to 197 (2015)
#' @keywords package
#' @examples
#' 
#'     if(interactive())
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
#'    }
#' 
NULL





#' example cleavage sites
#' 
#' An example data set containing cleavage sites (peaks) from getPeaks
#' 
#' 
#' @name peaks.gr
#' @docType data
#' @format \describe{GRanges with count (peak height), bg (local background),
#' SNratio (signal noise ratio), p-value, and option adjusted p-value }
#' @return \item{peaks.gr}{GRanges with count (peak height), bg (local
#' background), SNratio (signal noise ratio), p-value, and option adjusted
#' p-value }
#' @source http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1695644
#' @keywords datasets
#' @examples
#' 
#'     data(peaks.gr)
#'     names(peaks.gr)
#'     peaks.gr
#' 
NULL





#' example unique cleavage sites
#' 
#' An example data set containing cleavage sites with unique UMI, generated
#' from getUniqueCleavageEvents
#' 
#' 
#' @name uniqueCleavageEvents
#' @docType data
#' @return \describe{ \item{cleavage.gr }{Cleavage sites with one site per UMI
#' as GRanges with metadata column total set to 1 for each range}
#' \item{unique.umi.plus.R2}{a data frame containing unique cleavage site from
#' R2 reads mapped to plus strand with the following columns chr.y (chromosome
#' of readSide.y/R2 read) chr.x (chromosome of readSide.x/R1 read) strand.y
#' (strand of readSide.y/R2 read) strand.x (strand of readSide.x/R1 read)
#' start.y (start of readSide.y/R2 read) end.x (start of readSide.x/R1 read)
#' UMI (unique molecular identifier (umi) or umi with the first few bases of R1
#' read) } \item{unique.umi.minus.R2}{a data frame containing unique cleavage
#' site from R2 reads mapped to minus strand with the following columns chr.y
#' (chromosome of readSide.y/R2 read) chr.x (chromosome of readSide.x/R1 read)
#' strand.y (strand of readSide.y/R2 read) strand.x (strand of readSide.x/R1
#' read) end.y (end of readSide.y/R2 read) start.x (start of readSide.x/R1
#' read) UMI (unique molecular identifier (umi) or umi with the first few bases
#' of R1 read) } \item{unique.umi.plus.R1}{a data frame containing unique
#' cleavage site from R1 reads mapped to minus strand without corresponding R2
#' reads mapped to the plus strand, with the following columns chr.y
#' (chromosome of readSide.y/R2 read) chr.x (chromosome of readSide.x/R1 read)
#' strand.y (strand of readSide.y/R2 read) strand.x (strand of readSide.x/R1
#' read) start.x (start of readSide.x/R1 read) start.y (start of readSide.y/R2
#' read) UMI (unique molecular identifier (umi) or umi with the first few bases
#' of R1 read) } \item{unique.umi.minus.R1}{a data frame containing unique
#' cleavage site from R1 reads mapped to plus strand without corresponding R2
#' reads mapped to the minus strand, with the following columns chr.y
#' (chromosome of readSide.y/R2 read) chr.x (chromosome of readSide.x/R1 read)
#' strand.y (strand of readSide.y/R2 read) strand.x (strand of readSide.x/R1
#' read) end.x (end of readSide.x/R1 read) end.y (end of readSide.y/R2 read)
#' UMI (unique molecular identifier (umi) or umi with the first few bases of R1
#' read) } \item{all.umi}{a data frame containing all the mapped reads with the
#' following columns.  readName (read ID), chr.x (chromosome of readSide.x/R1
#' read), start.x (start of eadSide.x/R1 read), end.x (end of eadSide.x/R1
#' read), mapping.qual.x (mapping quality of readSide.x/R1 read), strand.x
#' (strand of readSide.x/R1 read), cigar.x (CIGAR of readSide.x/R1 read) ,
#' readSide.x (1/R1), chr.y (chromosome of readSide.y/R2 read) start.y (start
#' of readSide.y/R2 read), end.y (end of readSide.y/R2 read), mapping.qual.y
#' (mapping quality of readSide.y/R2 read), strand.y (strand of readSide.y/R2
#' read), cigar.y (CIGAR of readSide.y/R2 read), readSide.y (2/R2) R1.base.kept
#' (retained R1 length), R2.base.kept (retained R2 length), distance (distance
#' between mapped R1 and R2), UMI (unique molecular identifier (umi) or umi
#' with the first few bases of R1 read) }}
#' @source http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1695644
#' @keywords datasets
#' @examples
#' 
#'     data(uniqueCleavageEvents)
#'     names(uniqueCleavageEvents)
#'     sapply(uniqueCleavageEvents, class)
#'     uniqueCleavageEvents[[1]]  # GRanges object
#'     lapply(uniqueCleavageEvents, dim)
#' 
NULL




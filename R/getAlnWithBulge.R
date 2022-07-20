#' @importFrom Biostrings DNAStringSet nucleotideSubstitutionMatrix pairwiseAlignment readDNAStringSet reverseComplement
#' @importFrom utils read.table
#' @author Lihua Julie Zhu

getAlnWithBulge <- function(gRNA, gRNA.name,
                            peaks,
                            BSgenomeName,
                            mat,
                            peaks.withHeader = FALSE,
                            peaks.format= "bed",
                            gapOpening = 1L,
                            gapExtension = 3L,
                            max.DNA.bulge = 2L,
                            max.mismatch = 10L,
                            allowed.mismatch.PAM = 2L,
                            upstream = 20L,
                            downstream =20L,
                            PAM.size = 3L,
                            gRNA.size = 20L,
                            PAM = "NGG",
                            PAM.pattern = "NNN$",
                            PAM.location = "3prime")
{
  if (missing(gRNA.name))
    stop("gRNA name is required, e.g., testgRNA!")
  if (missing(gRNA) || nchar(gRNA) != gRNA.size)
    stop("gRNA is a required input with length equal to the specified
         gRNA.size, e.g., GACCCCCTCCACCCCGCCTC")
  if (missing(mat))
  {
    mat <- nucleotideSubstitutionMatrix(match = 1,
                                        mismatch = -1,
                                        baseOnly = TRUE)
    mat <- cbind(mat, N = -6)
    mat <- rbind(mat, N = -6)
  }
  if (missing(peaks))
    stop("peak is either DNAStringSet object or a path to a file in bed or
         fasta file format")
  if (class(peaks) != "DNAStringSet")
  {
    if (peaks.format == "bed")
    {
        if (missing(BSgenomeName))
        {
              stop("BSgenomeName is required if peaks is in bed file format!")
        }
        subjects2 <- CRISPRseek:::getSeqFromBed(inputFilePath = peaks,
                                              header = peaks.withHeader,
                                              BSgenomeName = BSgenomeName,
                                              upstream = upstream,
                                           downstream = downstream)
    }
    else if (peaks.format == "fasta")
    {
      subjects2 <- readDNAStringSet(peaks, format = format,
                                   use.names = TRUE)
    }
    else
    {
      stop("Only support DNAStringSet object, or files in bed, fasta, or fastq format.")
    }
  }
  else
  {
    subjects2 <- peaks
  }
  outfile <- tempfile(tmpdir = getwd())
  seqname <- names(subjects2)
  seqname <- gsub("'", "", seqname)
  seqname <- gsub(" ", "", seqname)
  seqname <- gsub("\t", ":", seqname)
  names(subjects2) <- seqname

  PAM.pattern <- gsub("^", "", gsub("$", "", PAM.pattern, fixed = TRUE),
                      fixed = TRUE)

  if (PAM.pattern != paste(rep("N", PAM.size), collapse=""))
  {
    #need to obtain separate masked.seq for forward and reverse strands
    stop("PAM pattern other than NNN is not fully implemented yet.")
    subSeqs <- .PAMpatternSearch(PAM.pattern = PAM.pattern,
                                    sequences = as.character(subjects2),
                                    PAM.location = PAM.location,
                                    PAM.size = PAM.size,
                                    gRNA.size = gRNA.size,
                                    max.DNA.bulge = max.DNA.bulge)
    seq.f <- subset(subSeqs,
              subSeqs$offTargetStrand == "+")

    seq.r <- subset(subSeqs,
              subSeqs$offTargetStrand == "-")

    seq.f.DSS <- DNAStringSet(seq.f$ProtoSpacer)
    seq.r.DSS <- DNAStringSet(seq.r$ProtoSpacer)

    pa.f.1 <- lapply(1:length(seq.f.DSS), function(i) {
      pairwiseAlignment(pattern = gRNA,
                        subject = seq.f.DSS[i],
                        type = "global-local",
                        substitutionMatrix = mat,
                        gapOpening = gapOpening,
                        gapExtension = gapExtension,
                        scoreOnly = FALSE)
    })

    # pa.r.2 <- lapply(1:length(seq.r.DSS), function(i) {
    #   pairwiseAlignment(pattern =
    #                       as.character(reverseComplement(DNAString(gRNA))),
    #                     subject = seq.r.DSS[i],
    #                     type = "global-local",
    #                     substitutionMatrix = mat,
    #                     gapOpening = gapOpening,
    #                     gapExtension = gapExtension,
    #                     scoreOnly = FALSE)
    # })

    pa.r.2 <- lapply(1:length(seq.r.DSS), function(i) {
      pairwiseAlignment(pattern =
                         gRNA,
                        subject = as.character(reverseComplement(seq.r.DSS[i])),
                        type = "global-local",
                        substitutionMatrix = mat,
                        gapOpening = gapOpening,
                        gapExtension = gapExtension,
                        scoreOnly = FALSE)
    })
    # need to change subjects2 to multiple of fragments if choosing multiple
    # fragments instead of masking as an alternative strategy
    best.aln.info <- do.call(rbind, lapply(1:length(seq.f.DSS), function(i) {
      getBestAlnInfo(seq.f.DSS[[i]], pa.f1[[i]], NA, PAM.size = PAM.size,
                     PAM = PAM, gRNA.size = gRNA.size)
    }))

    best.aln.info <- do.call(rbind, lapply(1:length(seq.r.DSS), function(i) {
      getBestAlnInfo(seq.r.DSS[i], NA, pa.r2[[i]], PAM.size = PAM.size,
                     PAM = PAM, gRNA.size = gRNA.size)
    }))
  }
  else
  {
    pa.f<- lapply(1:length(seqname), function(i) {
      pairwiseAlignment(pattern = gRNA,
                        subject = as.character(subjects2[i]),
                        type = "global-local",
                        substitutionMatrix = mat,
                        gapOpening = gapOpening,
                        gapExtension = gapExtension,
                        scoreOnly = FALSE)
    })
    # pa.r <- lapply(1:length(seqname), function(i) {
    #   pairwiseAlignment(pattern = as.character(
    #     reverseComplement(DNAString(gRNA))),
    #     subject = as.character(subjects2[i]),
    #     type = "global-local",
    #     substitutionMatrix = mat,
    #     gapOpening = gapOpening,
    #     gapExtension = gapExtension,
    #     scoreOnly = FALSE)
    #})
     # favor matches with bulge distant from PAM
    pa.r <- lapply(1:length(seqname), function(i) {
        pairwiseAlignment(pattern = gRNA,
          subject = as.character(reverseComplement(subjects2[i])),
          type = "global-local",
          substitutionMatrix = mat,
          gapOpening = gapOpening,
          gapExtension = gapExtension,
          scoreOnly = FALSE)
      })

    best.aln.info <- do.call(rbind, lapply(1:length(seqname), function(i) {
      getBestAlnInfo(subjects2[i], pa.f[[i]], pa.r[[i]], PAM.size = PAM.size,
                     PAM = PAM, gRNA.size = gRNA.size)
    }))
  }
  best.aln.info <- as.data.frame(best.aln.info)
  best.aln.info$gRNA.name = gRNA.name
  if (PAM.location == "3prime")
    best.aln.info$gRNAPlusPAM = paste0(gRNA, PAM)
  else
    best.aln.info$gRNAPlusPAM = paste0(PAM, gRNA)

  if (class(peaks) != "DNAStringSet" && peaks.format == "bed")
  {
    bed <- read.table(peaks, sep = "\t", header = peaks.withHeader,
                      stringsAsFactors = FALSE)
    if (dim(bed)[2] >= 5)
      colnames(bed)[5] = "peak_score"
    else
      bed$peak_score <- NA

    if (dim(bed)[2] >= 6) {
      strand <- bed[, 6]
      strand[!strand %in% c("+", "-", "*")] <- "*"
    }
    else {
      strand <- rep("+", dim(bed)[1])
    }
    chr <- bed[, 1]
    Start <- as.numeric(bed[, 2])
    End <- as.numeric(bed[, 3])
    Start <- ifelse(strand == "-", Start - downstream, Start -
                      upstream)
    End <- ifelse(strand == "-", End + upstream, End + downstream)
    offTarget_Start <- unlist(apply(cbind(Start, 1), 1, max))

    best.aln.info$offTarget_Start <- unlist(best.aln.info$offTarget_Start) +
                                        offTarget_Start - 1
    best.aln.info$offTarget_End <- unlist(best.aln.info$offTarget_End) +
                                        offTarget_Start - 1
    best.aln.info$chromosome <- chr
    best.aln.info$peak_score <- bed$peak_score
    best.aln.info$offTarget <- bed[,4]
  }

  n.indel <- unlist(best.aln.info$n.insertion) +
    unlist(best.aln.info$n.deletion)

  list(aln.all = best.aln.info,
       aln.indel = best.aln.info[n.indel > 0 &
                    best.aln.info$n.PAM.mismatch <= allowed.mismatch.PAM &
                    (n.indel +
                       unlist(best.aln.info$n.mismatch)) <=
                    max.mismatch,])
}

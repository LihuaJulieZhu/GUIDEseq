#' Parse pairwise alignment
#'
#' @param offtargetSeq DNAStringSet object of length 1
#' @param pa.f Global-Local PairwiseAlignmentsSingleSubject, results of
#' pairwiseAlignment, alignment of pattern to subject
#' @param pa.r Global-Local PairwiseAlignmentsSingleSubject, results of
#' pairwiseAlignment, alignment of reverse pattern to subject
#' @param gRNA.size size of gRNA, default to 20
#' @param PAM PAM sequence, default to NGG
#' @param PAM.size PAM size, default to 3
#' @param insertion.symbol symbol for representing bulge in offtarget,
#' default to ^. It can also be set to lowerCase to use lower case letter
#' to represent insertion
#'
#' @return a dataframe with the following columns.
#' offTarget: name of the offtarget
#' peak_score: place holder for storing peak score
#' predicted_cleavage_score: place holder for storing cleavage score
#' gRNA.name: place holder for storing gRNA name
#' gRNAPlusPAM: place holder for storing gRNAPlusPAM sequence
#' offTarget_sequence: offTarget sequence with PAM in the right orientation.
#' For PAM in the 3' prime location, offTarget is the sequence on the plus strand
#' otherwise, is the sequence on the reverse strand
#' seq.aligned: the aligned sequence without PAM
#' guideAlignment2OffTarget: string representation of the alignment
#' offTargetStrand: the strand of the offtarget
#' mismatch.distance2PAM: mismatch distance to PAM start
#' n.PAM.mismatch: number of mismatches in PAM
#' n.guide.mismatch: number of mismatches in the gRNA not including PAM
#' PAM.sequence: PAM in the offtarget
#' offTarget_Start: offtarget start
#' offTarget_End: offTarget end
#' chromosome: place holder for storing offtarget chromosome
#' pos.mismatch: mismatch positions with the correct PAM orientation, i.e.,
#' indexed form distal to proximal of PAM
#' pos.indel: indel positions starting with deletions in the gRNA followed
#' by those in the offtarget
#' pos.insertion: Insertion positions in the gRNA
#' Insertion positions are counted from distal to proximal of PAM
#' For example, 5 means the 5th position is an insertion in
#' gRNA
#' pos.deletion: Deletion in the gRNA
#' Deletion positions are counted from distal to proximal of PAM
#' For example, 5 means the 5th position is a deletion in
#' gRNA
#' n.insertion: Number of insertions in the RNA. Insertions in gRNA creates
#' bulged DNA base
#' n.deletion: Number of deletions in the RNA. Deletions in gRNA creates
#' bulged DNA base
#' @importFrom Biostrings score start alignedPattern alignedSubject reverseComplement
#' DNAString subject
#' @importFrom Base names substr which class missing is.na strsplit length is.na
#' @export
#'
#' @examples
getBestAlnInfo <- function(offtargetSeq, pa.f, pa.r, gRNA.size = 20,
                           PAM = "NGG", PAM.size = 3, insertion.symbol = "^")
{
  if (missing(offtargetSeq))
    stop("offtargetSeq is required!")
  if (is.na(pa.f) && is.na(pa.r))
    stop("At least one of pairwise alignment objects pa.f or pa.r is not NA!")

  if (!is.na(pa.f) && !is.na(pa.r) &&
      class(pa.f) == "PairwiseAlignmentsSingleSubject" &&
      class(pa.r) == "PairwiseAlignmentsSingleSubject")
  {
    scores <- c(score(pa.f),
                score(pa.r))
    match.starts <- c(start(subject(pa.f)),
                      start(subject(pa.r)))
    best.start <- match.starts[scores == max(scores)][1]
    best.pos <- which(scores == max(scores))[1]
    best.aln <- switch(best.pos,
                       pa.f,
                       pa.r
    )
  }
  else if (!is.na(pa.f) && class(pa.f) == "PairwiseAlignmentsSingleSubject")
  {
    best.aln <- pa.f
    scores <- score(pa.f)
    match.starts <- start(subject(pa.f))
    best.start <- match.starts[scores == max(scores)][1]
    best.pos <- which(scores == max(scores))[1]
  }
  else if (!is.na(pa.r) && class(pa.r) == "PairwiseAlignmentsSingleSubject")
  {
    best.aln <- pa.r
    scores <- score(pa.r)
    match.starts <- start(subject(pa.r))
    best.start <- match.starts[scores == max(scores)][1]
    best.pos <- which(scores == max(scores))[1]
  }
 else
 {
   stop("pa.f or pa.r is a required input as PairwiseAlignmentsSingleSubject")
 }

  seq <- c(alignedPattern(best.aln), alignedSubject(best.aln))

  name.peak <- names(seq)[2]

  seq.print <- as.character(
    switch(best.pos[1],
           seq[2],
           reverseComplement(DNAString(as.character(seq[2])))))

  pos.mismatch <- switch(best.pos,
                         pa.f@subject@mismatch[[1]] - best.start + 1,
                         nchar(seq[1]) - (pa.r@subject@mismatch[[1]] -
                                        best.start[1]))

  guide.print <- as.character(seq[1])
  seq.print<- strsplit(seq.print, "")[[1]]

  guide.print.v <- strsplit(guide.print, "")[[1]]
  if(best.pos == 2)
  {
      pos.deletion.guide <- gRNA.size - (which(guide.print.v == "-") - 1) + 1
   }
  else
  {
      pos.deletion.guide <- which(guide.print.v == "-")
  }
  pos.deletion.off <- which(seq.print == "-")
  if(insertion.symbol == "lowerCase")
    seq.print[pos.deletion.guide] <- tolower(seq.print[pos.deletion.guide])
  else
    seq.print[pos.deletion.guide] <- insertion.symbol

  if (length(pos.deletion.off) > 0  )
  {
    for (i in 1:length(pos.deletion.off))
    {
      if (best.pos == 1)
        pos.mismatch[(pos.mismatch - pos.deletion.off[i]) > 0] <-
          pos.mismatch[(pos.mismatch - pos.deletion.off[i]) > 0] + 1
      else
        pos.mismatch[(pos.mismatch - pos.deletion.off[i]) < 0] <-
           pos.mismatch[(pos.mismatch - pos.deletion.off[i]) < 0] - 1
    }
  }
 # It is important to set seq.print before reset pos.mismatch using
 # pos.deletion.guide information

  n.deletion <- length(which(guide.print.v == "-"))
  n.insertion <- length(which(seq.print == "-"))

  seq.print[setdiff(1:width(seq[2]),c(pos.deletion.guide,
                                      pos.deletion.off,
                                      pos.mismatch))] <- "."


  if (length(pos.deletion.guide) > 0  )
  {
    for (i in 1:length(pos.deletion.guide))
    {
      if (best.pos == 1)
        pos.mismatch[(pos.mismatch - pos.deletion.guide[i]) > 0] <-
          pos.mismatch[(pos.mismatch - pos.deletion.guide[i]) > 0] - 1
      else
      {
        pos.mismatch[(pos.mismatch - pos.deletion.guide[i]) < 0] <-
          pos.mismatch[(pos.mismatch - pos.deletion.guide[i]) < 0] + 1
        # need to add the following line after seq.print is set
        pos.mismatch <- pos.mismatch - nchar(seq[1]) + gRNA.size
      }
    }
  }

  seq.aligned <- switch(best.pos,
                        as.character(pa.f@subject),
                        as.character(reverseComplement(DNAString(as.character(
                          pa.r@subject)))))
  strand.aligned <- switch(best.pos,
                           "+", "-")

  PAM.aligned <- switch(best.pos,
                        substr(offtargetSeq, best.start + width(seq[2]) -
                                 length(which(seq.print == "-")) ,
                               best.start + width(seq[2]) + PAM.size - 1 -
                                 length(which(seq.print == "-"))),
                        as.character(reverseComplement(DNAString(
                          substr(offtargetSeq, best.start - PAM.size,
                                 best.start - 1))))
  )
  seq.aligned <- paste0(seq.aligned, PAM.aligned)

  # only allow A, C, G, T, and N in the PAM sequence for now.
  # need to add other code types later

  offTarget_Start <- switch(best.pos,
                           best.start,
                           max(1, best.start - PAM.size))

  offTarget_End <- offTarget_Start + width(seq[2]) - 1 + PAM.size -
    n.insertion

 seq.print<- paste0(seq.print, collapse = "")

 list(offTarget = name.peak,
      peak_score = NA,
      predicted_cleavage_score = NA,
      gRNA.name = NA,
      gRNAPlusPAM = NA,
      offTarget_sequence = seq.aligned,
      guideAlignment2OffTarget = seq.print,
      offTargetStrand = strand.aligned,
      mismatch.distance2PAM = paste(gRNA.size - pos.mismatch + 1,
                                    collapse = ","),
      n.PAM.mismatch	=  neditAt(DNAString(PAM.aligned),
                                DNAString(PAM), fixed=FALSE),
      n.mismatch = length(pos.mismatch),
      PAM.sequence = PAM.aligned,
      offTarget_Start = offTarget_Start,
      offTarget_End = offTarget_End,
      chromosome = NA,
      pos.mismatch = pos.mismatch,
      pos.insertion = pos.deletion.off,
      pos.deletion = pos.deletion.guide,
      n.insertion = n.insertion,
      n.deletion = n.deletion)
      #seq.aln = seq)
}

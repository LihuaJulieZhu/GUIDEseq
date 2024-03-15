#' @importFrom CRISPRseek translatePattern

.addPAMpattern <- function(gRNA, PAM.pattern, PAM.size, PAM.location)
{
  PAM.pattern <- gsub("^", "", gsub("$", "", PAM.pattern, fixed = TRUE),
                      fixed = TRUE)
  if (PAM.pattern != paste(rep("N", PAM.size), collapse=""))
  {
    if ( PAM.location == "3prime" && nchar(PAM.pattern) == PAM.size)
      gRNAs <- paste0(gRNA, PAM.pattern)
    else if  (PAM.location == "5prime" && nchar(PAM.pattern) == PAM.size)
      gRNAs <- paste0(PAM.pattern, gRNA)
    else if (length(strsplit(PAM.pattern, "|")[[1]]) > 1)
    {
      PAM.patterns <- strsplit(PAM.pattern, "|", fixed = TRUE)[[1]]
      gRNAs <- unlist(lapply(PAM.patterns, function(this.pattern) {
        .addPAMpattern(gRNA, this.pattern, PAM.size, PAM.location)
      }))
    }
    else
      stop("Please specify PAM.pattern using IUPAC nucleotide ambiguity codes!")
  }
  else
    gRNAs <- gRNA
  gRNAs
}


.getPAMInd <- function(PAM.pattern, sequences)
{
  PAM.pattern <- gsub("^", "",
                      gsub("$", "", PAM.pattern, fixed = TRUE),
                      fixed = TRUE)
  if (length(strsplit(PAM.pattern, "|", fixed = TRUE)[[1]]) <= 1)
    PAM.pattern <- translatePattern(PAM.pattern)
  PAM.pattern <- paste0("(?=", PAM.pattern, ")")
  pattern.ind <- gregexpr(PAM.pattern, sequences,
                          perl = TRUE, ignore.case = TRUE)
  pattern.ind
}

.getMatchedPAMInd <- function(inds,
                           direction = c("forward", "reverse"),
                           sequences, PAM.size = 3L,
                           gRNA.size = 20, PAM.location = "3prime")
{
  inds <- as.numeric(unlist(inds))
  direction <- match.arg(direction)

  if (PAM.location == "3prime")
  {
    if (direction == "forward")
    {
      ind.start <- inds - gRNA.size
      ind.end <- inds - 1
      ind.end[ind.end < 1] <- 1
      ind.start[ind.start < 1] <- 1
    }
    else if(max(as.numeric(inds)) > 0)
    {
      ind.start <- inds + PAM.size
      ind.end <- inds + PAM.size + gRNA.size - 1
      ind.end[ind.end > nchar(sequences)] <- nchar(sequences)
      ind.end[ind.end < 1] <- 1
    }
    else
    {
      ind.start <- 1
      ind.end <- 1
    }
    if (max(as.numeric(inds)) > 0 &&
        length(which((ind.end - ind.start + 1) < gRNA.size )))
      warning("PAM pattern found near the edge of the search region.
              You may need to increase upstream and downstream to include
              a wider search region!")
     list(ind.start = ind.start, ind.end = ind.end)
  }
  else
  {
    stop("Sorry! 5prime has not been implemented yet!")
  }
}

.getSubSeq <- function(inds, this.seq, direction = c("forward", "reverse"),
                       PAM.size = 3L,
                       gRNA.size = 20L,
                       PAM.location = "3prime",
                       max.DNA.bulge = 2L)
{
  if (PAM.location == "3prime")
  {
    temp <- .getMatchedPAMInd(inds, this.seq, direction = direction,
                           PAM.size = PAM.size,
                           gRNA.size = gRNA.size,
                           PAM.location = PAM.location)

    ind.start <- temp[[1]] # start of prospacer not including PAM
    ind.end <- temp[[2]] # end of prospacer (not including PAM)
    if (direction == "forward")
             sub.seq <- do.call(rbind, lapply(1:length(ind.start),
                               function(i) {
                                c(this.seq, "+", max(1, ind.start[i] - max.DNA.bulge),
                                  ind.end[i] + PAM.size,
                                  ifelse (ind.end[i] == 1, NA,
                                            substr(this.seq,
                                                  max(1, ind.start[i] - max.DNA.bulge),
                                            ind.end[i])),
                                  ifelse (ind.end[i] == 1, NA,
                                            substr(this.seq, ind.end[i] + 1,
                                                       ind.end[i] + PAM.size)))}))
    else
            sub.seq <- do.call(rbind, lapply(1:length(ind.start),
                                function(i) {
                                     c(this.seq, "-", ind.start[i] - PAM.size,
                                       min(nchar(this.seq),
                                           ind.end[i] + max.DNA.bulge),
                                       ifelse (ind.end[i] == 1, NA,
                                               as.character(reverseComplement(
                                                 DNAString(substr(this.seq,
                                                     ind.start[i],
                                                      min(nchar(this.seq),
                                                          ind.end[i] + max.DNA.bulge)))))),
                                      ifelse (ind.end[i] == 1, NA,
                                               as.character(reverseComplement(
                                                 DNAString(substr(this.seq,
                                                       ind.start[i] - PAM.size,
                                                       ind.start[i] - 1))))))}))
    sub.seq <- subset(sub.seq, !is.na(sub.seq[,5]))
    temp <- cbind(sub.seq, paste0(sub.seq[,5], sub.seq[,6]))
    colnames(temp) <- c("Seq2Search", "offTargetStrand",
                        "offTarget_Start", "offTarget_End",
                        "ProtoSpacer", "PAM", "offTarget_sequence")
    temp

  }
  else
  {
    stop("5prime has not been implemented yet!")
  }
}

.PAMpatternSearch <- function(PAM.pattern, sequences,
                              PAM.location = "3prime",
                              PAM.size = 3,
                              gRNA.size = 20,
                              max.DNA.bulge = 2L)
{
  # For minus strand, offTarget_sequence is the reverse complement of the input
  # sequence with PAM in the right orientation
  # For 3prime PAM, offTarget_Start and offTarget_Ends are
  # on the plus strand: PAM starts at offTarget_Ends - PAM.size + 1 and ends at offTarget_Ends
  # on the minus strand: PAM starts at PAM.size and ends at 1
  # Seq2Search is the input sequences, each input sequence could have 0 to many
  # potential offtarget sites with the PAM.pattern

  PAM.pattern <- gsub("^", "",
                      gsub("$", "", PAM.pattern, fixed = TRUE),
                      fixed = TRUE)
  rev.pattern <- as.character(reverseComplement(DNAString(PAM.pattern)))

  ind.f <- .getPAMInd(PAM.pattern = PAM.pattern, sequences = sequences)
  ind.r <- .getPAMInd(PAM.pattern = rev.pattern , sequences = sequences)

  sub.sequences.f <- do.call(rbind, lapply(1:length(sequences), function(i) {
     .getSubSeq(ind.f[[i]],
               sequences[i], direction = "forward",
               PAM.size = PAM.size,
               gRNA.size = gRNA.size,
               PAM.location  = PAM.location,
               max.DNA.bulge = max.DNA.bulge)
   }))

   sub.sequences.r <- do.call(rbind, lapply(1:length(sequences), function(i) {
     .getSubSeq(ind.r[[i]],
                sequences[i], direction = "reverse",
                PAM.size = PAM.size,
                gRNA.size = gRNA.size,
                PAM.location  = PAM.location,
                max.DNA.bulge = max.DNA.bulge)

   }))
   as.data.frame(rbind(sub.sequences.f, sub.sequences.r))
  # unlist(lapply(1:length(sequences), function(i) {
  #   .maskSubSeq(ind.f[[i]], ind.r[[i]], sequences[[i]],
  #               PAM.location = PAM.location)
  # }))
}


.maskSubSeq <- function(ind.fi, ind.ri, .sequence, PAM.size = 3L,
                        gRNA.size = 20L,
                        PAM.location = "3prime")
{
  if (PAM.location == "3prime")
  {
    temp <- .getMatchedInd(ind.fi, ind.ri, .sequence,
                           PAM.size = PAM.size,
                           gRNA.size = gRNA.size,
                           PAM.location = PAM.location)

    .getMaskedSeq <- function(ind.start, ind.end)
    {
      gap.range <- gaps(reduce(IRanges(ind.start, ind.end)),
                        start = 1, end = nchar(.sequence))
      temp.seq <- .sequence
      for (i in 1:length(gap.range))
      {
        substr(temp.seq, start(gap.range)[i], end(gap.range)[i]) <-
          paste(rep("N", width(gap.range[i])), collapse="")
      }
      temp.seq
    }
    list(masked.f.seq = .getMaskedSeq(temp[[1]][[1]], temp[[1]][[2]]),
         masked.r.seq = .getMaskedSeq(temp[[2]][[1]], temp[[2]][[2]]))
  }
  else
  {
    stop("Sorry! 5prime not implemented yet!")
  }
}

.nucleotideSubstitutionMatrix <- function(match = 1,
                                          mismatch = 0,
                                          baseOnly = FALSE,
                                          type = "DNA")
{
  "%safemult%" <- function(x, y) ifelse(is.infinite(x) & y ==
                                          0, 0, x * y)
  type <- match.arg(type, c("DNA", "RNA"))
  if (!isSingleNumber(match) || !isSingleNumber(mismatch))
    stop("'match' and 'mismatch' must be non-missing numbers")
  if (baseOnly)
    letters <- IUPAC_CODE_MAP[DNA_BASES]
  else letters <- IUPAC_CODE_MAP
  if (type == "RNA")
    names(letters) <- chartr("T", "U", names(letters))
  nLetters <- length(letters)
  splitLetters <- strsplit(letters, split = "")
  submat <- matrix(0, nrow = nLetters, ncol = nLetters,
                   dimnames = list(names(letters),
                                   names(letters)))
  for (i in 1:nLetters)
  {
    for (j in i:nLetters)
    {
      if (i == j ||
          i %in% which(rownames(submat) %in% unlist(splitLetters[colnames(submat)[j]])))
        submat[i, j] <- submat[j, i] <- match
      else
        submat[i, j] <- submat[j,i] <- mismatch
    }
  }
  submat
  #mean(outer(splitLetters[[i]], splitLetters[[j]], "=="))
  #abs(match) * submat - abs(mismatch) %safemult% (1 - submat)
}

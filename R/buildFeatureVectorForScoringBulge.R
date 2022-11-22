.mismatches_as_IntegerList <- function(mismatches)
{
    arr_ind <- which(mismatches != 0, arr.ind=TRUE)
    ind_row <- unname(arr_ind[ , "row"])
    ind_col <- unname(arr_ind[ , "col"])
    oo <- S4Vectors::orderIntegerPairs(ind_row, ind_col)
    ind_row <- ind_row[oo]
    ind_col <- ind_col[oo]
    partitioning <- PartitioningByEnd(ind_row, NG=nrow(mismatches))
    relist(ind_col, partitioning)
}

#' Build Feature Vector For Scoring Offtargets with Bulge
#'
#' @param alns alignments, output from getAlnWithBulge
#' (see the example below)
#'
#' @param gRNA.size Size of the gRNA, default to 20L
#' @param canonical.PAM PAM sequence, default to NGG
#' @param subPAM.start  start of the subPAM, default to 2L for NGG
#' @param subPAM.end    End of the subPAM, default to 3L for NGG
#' @param insertion.symbol Symbol used to indicate bulge in DNA
#' Default to ^
#' @param PAM.size Size of the PAM, default to 3L for NGG
#' @param PAM.location The location of the PAM, default to 3prime
#'
#' @author Lihua Julie Zhu

#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' isMatchingAt substring extractAt isMatchingAt
#' complement DNAString replaceAt
#' @importFrom BiocGenerics grep sort relist rep.int cbind
#' @importFrom IRanges IRangesList PartitioningByEnd
#' @importFrom S4Vectors unstrsplit orderIntegerPairs elementNROWS
#' @importFrom methods as
#' @importFrom BiocGenerics relist
#'
#' @export buildFeatureVectorForScoringBulge
#'
#' @examples
#'
#' if (interactive())
#' {
#'   library(BSgenome.Hsapiens.UCSC.hg19)
#'   library(GUIDEseq)
#'   peaks.f <- system.file("extdata", "T2plus100OffTargets.bed",
#'      package = "GUIDEseq")
#'   gRNA <- "GACCCCCTCCACCCCGCCTC"
#'   temp <- GUIDEseq:::getAlnWithBulge(gRNA, gRNA.name = "T2",
#'       peaks = peaks.f, BSgenomeName = Hsapiens,
#'        peaks.withHeader = TRUE)
#'    fv <- buildFeatureVectorForScoringBulge(temp$aln.indel)
#'    fv$featureVectors
#'  }

buildFeatureVectorForScoringBulge <-
    function(alns, gRNA.size = 20,
    canonical.PAM = "NGG", subPAM.start = 2, subPAM.end = 3,
    insertion.symbol = "^",
     PAM.size = 3, PAM.location = "3prime")
{
    if (PAM.location != "3prime")
        stop("Bulge scoring for 5prime PAM not implemented yet")
    if (dim(alns)[1] == 0)
    {
        stop("Empty hits!")
    }
    #subject <- DNAStringSet(as.character(alns$offTarget_sequence))

    # remove alignments with extended peak sequences too short
    # users should increase upstream and downstream to fetch sufficient number
    # of bases for givent peaks if needed
      
    alns <- alns[alns$PAM != "" & nchar(alns$PAM) == PAM.size, ]
    
    PAM <- DNAStringSet(unlist(alns$PAM.sequence))
    isCanonical.PAM <- as.numeric(isMatchingAt(canonical.PAM,
                                               PAM,
                                               at = 1, fixed = FALSE))
    mismatch_pos <- alns$pos.mismatch
    insertion_pos <- alns$pos.insertion
    deletion_pos <- alns$pos.deletion

    mismatches <- matrix(rep(0, gRNA.size * length(PAM)),
                         ncol = gRNA.size, nrow = length(PAM))
    colnames(mismatches) <- paste0("IsMismatch.pos", 1:gRNA.size)
    for (i in 1:nrow(mismatches)) {
        if (length(mismatch_pos[[i]][mismatch_pos[[i]] <=
                gRNA.size & mismatch_pos[[i]] > 0]) > 0)
          mismatches[i, mismatch_pos[[i]][mismatch_pos[[i]] <= 
                gRNA.size & mismatch_pos[[i]] > 0]] <- 1
    }
    insertions <- matrix(rep(0, gRNA.size * length(PAM)),
                         ncol = gRNA.size, nrow = length(PAM))
    colnames(insertions) <- paste0("IsInsertion.pos", 1:gRNA.size)
    for (i in 1:nrow(insertions)) {
      if (length(insertion_pos[[i]][insertion_pos[[i]] <=
                gRNA.size & insertion_pos[[i]] > 0]) >0 )
          insertions[i, insertion_pos[[i]][insertion_pos[[i]] <=
                gRNA.size & insertion_pos[[i]] > 0]] <- 1
    }
    deletions <- matrix(rep(0, gRNA.size * length(PAM)),
                         ncol = gRNA.size, nrow = length(PAM))
    colnames(deletions) <- paste0("IsDeletion.pos", 1:gRNA.size)
    for (i in 1:nrow(deletions)) {
      if (length(deletion_pos[[i]][deletion_pos[[i]] <=
                gRNA.size & deletion_pos[[i]] > 0]) > 0)
          deletions[i, deletion_pos[[i]][deletion_pos[[i]] <=
                gRNA.size & deletion_pos[[i]] > 0]] <- 1
    }
    #d.nucleotide is the nucleotide that will hybridize to the gRNA
    ### reverse complement of the offtarget sequence
    #r.nucleotide is the gRNA sequence,except T is converted to U)

    at <- IRangesList(start=mismatch_pos, end=mismatch_pos)
    at.ins <- IRangesList(start=insertion_pos, end=insertion_pos)
    at.del <- IRangesList(start=deletion_pos, end=deletion_pos)

    if (insertion.symbol == "lowerCase")
    {
        insertion.symbol =="[acgt]"
        OT <- gsub(".", "A", gsub(insertion.symbol, "",
                   alns$guideAlignment2OffTarget,
                   fixed = FALSE), fixed = TRUE)
    }
    else
    {
        OT <- gsub(".", "A", gsub(insertion.symbol, "",
                   alns$guideAlignment2OffTarget,
                   fixed = TRUE), fixed = TRUE)
    }

    r.ins <- as.character(unlist(extractAt(
      DNAStringSet(unlist(alns$gRNAPlusPAM)), at.ins)))
    r.del <- as.character(unlist(extractAt(
      DNAStringSet(unlist(alns$offTarget_sequence)), at.del)))
    # insertion and deletion are on gRNA, so need to translate from T to U
    # for CFD score calculation
    # need to convert the deletion back to T in the offtarget output file
    r.ins[r.ins == "T"] <- "U"
    r.del[r.del == "T"] <- "U"

    d.nucleotide <- extractAt(complement(DNAStringSet(OT)), at)
    r.nucleotide <- as.character(unlist(extractAt(
        DNAStringSet(alns$gRNAPlusPAM), at)))
    r.nucleotide[r.nucleotide == "T"] <- "U"
    d.nu.r.nu <- paste("r", r.nucleotide, ":d",
      unlist(d.nucleotide), sep="")
    #### need to assign d.nu.r.nu to the right alignment
    arr_ind <- which(mismatches != 0, arr.ind=TRUE)
    ind_row <- sort(unname(arr_ind[ , "row"]))
    partitioning <- PartitioningByEnd(ind_row, NG=nrow(mismatches))
    d.nu.r.nu.2 <- relist(d.nu.r.nu, partitioning)
    d.nu.r.nu <- unstrsplit(as(d.nu.r.nu.2,
                   "CharacterList"), sep = ",")

    #### need to assign insertion and deletion to the right alignment
    arr_ind <- which(insertions != 0, arr.ind=TRUE)
    ind_row <- sort(unname(arr_ind[ , "row"]))
    partitioning <- PartitioningByEnd(ind_row, NG=nrow(insertions))
    r.ins.2 <- relist(r.ins, partitioning)
    r.ins <- unstrsplit(as(r.ins.2,
                               "CharacterList"), sep = ",")

    arr_ind <- which(deletions != 0, arr.ind=TRUE)
    ind_row <- sort(unname(arr_ind[ , "row"]))
    partitioning <- PartitioningByEnd(ind_row, NG=nrow(deletions))
    r.del.2 <- relist(r.del, partitioning)
    r.del <- unstrsplit(as(r.del.2,
                           "CharacterList"), sep = ",")

    mismatch.distance2PAM <-  unlist(alns$mismatch.distance2PAM)

    # alignment <- rep.int(DNAString("."), gRNA.size)
    # alignment <- rep.int(DNAStringSet(alignment), nrow(hits))
    # alignment <- as.character(replaceAt(alignment, at, extractAt(subject, at)))

    mean.neighbor.distance.mismatch <- unlist(lapply(1:nrow(mismatches),
                                              function(i) {
                                                  mean(abs(diff(mismatch_pos[[i]])))
                                              }))
    no_neighbor_idx <- elementNROWS(mismatch_pos) <= 1L
    mean.neighbor.distance.mismatch[no_neighbor_idx] <- gRNA.size

    subPAM <- substr(alns$PAM.sequence, subPAM.start, subPAM.end)
    features <- cbind(isCanonical.PAM,
        mean.neighbor.distance.mismatch, d.nu.r.nu, subPAM, r.ins, r.del)
    #colnames(alns)[colnames(alns) == "guideAlignment2OffTarget"] <- "alignment"

    colnames(features) <- c("NGG",
        "mean.neighbor.distance.mismatch", "mismatch.type", "subPAM",
        "gRNA.insertion",
        "gRNA.deletion")
    temp <- cbind(alns, mismatches, insertions, deletions, features)
    subPAM.width <- nchar(subPAM)
    if (length(which(subPAM.width < subPAM.end - subPAM.start + 1)) > 0 )
    {
       retained.temp <- temp[subPAM.width == (subPAM.end - subPAM.start + 1),]
       warning("Some of the offtarget sites got dropped
       due to that the PAM sequence is outside of the searching region!
       Please increase the searching window by increasing upstream
               and downstream! For details, see the second list item in the
               returned object from buildFeatureVectorForScoringBulge")
       dropped.temp <- temp[subPAM.width < (subPAM.end - subPAM.start + 1),]
       list(featureVectors = retained.temp, featureVectors.dropped = dropped.temp)
    }
    else
    {
      list(featureVectors = temp, featureVectors.dropped = NA)
    }
}

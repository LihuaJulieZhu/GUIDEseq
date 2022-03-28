annotateOffTargets <-
function (thePeaks, txdb, orgAnn) 
{
    thePeaks <- subset(thePeaks, !is.na(offTarget_Start) & offTarget_Start != "")
    peaks.RD <- GRanges(seqnames = Rle(thePeaks$chromosome), 
        ranges = IRanges(start = thePeaks$offTarget_Start, 
        end = thePeaks$offTarget_End, names = thePeaks$names))
    allExons <- as(exons(txdb, columns = "gene_id"), "GRanges")
    if (seqlevelsStyle(allExons) != seqlevelsStyle(peaks.RD)) {
        seqlevelsStyle(allExons) <-  seqlevelsStyle(peaks.RD)
    }
    allExons <- allExons[as.character(seqnames(allExons)) %in% 
        unique(as.character(seqnames(peaks.RD))), ]
    inExon <- overlapsAny(peaks.RD, allExons, minoverlap = 1L, 
        type = "any", ignore.strand = TRUE)
    inExon[inExon== FALSE] <- ""
#    allGenes <- genes(txdb, vals = NULL, columns = "gene_id", 
     suppressMessages(allGenes <- genes(txdb, columns = "gene_id",
        single.strand.genes.only =  TRUE))
    if (seqlevelsStyle(allGenes) != seqlevelsStyle(peaks.RD)) {
        seqlevelsStyle(allGenes) <-  seqlevelsStyle(peaks.RD)
    }
    overlapGenes <- findOverlaps(peaks.RD, allGenes, minoverlap = 1L, 
        type = "any", ignore.strand = TRUE)
    entrez_id <- character(dim(thePeaks)[1])
    query.ind <- queryHits(overlapGenes)
    entrez_id[query.ind] <- unlist(allGenes[subjectHits(overlapGenes), ]$gene_id)
    entrez_id <- as.character(entrez_id)
    entrez_id[is.na(entrez_id)] = ""
    thePeaks <- cbind(thePeaks, inExon = inExon, entrez_id = entrez_id)
    if (length(queryHits(overlapGenes)) > 0 && !missing(orgAnn) && 
        class(orgAnn) == "AnnDbBimap") {
        egSYMBOL <- toTable(orgAnn)
        if (length(grep("flybase_id", colnames(egSYMBOL)[2])) > 
            0) {
            m <- match(thePeaks$entrez_id, egSYMBOL$flybase_id)
            thePeaks$symbol <- egSYMBOL[, 1][m]
        }
        else {
            m <- match(thePeaks$entrez_id, egSYMBOL$gene_id)
            thePeaks$symbol <- egSYMBOL$symbol[m]
        }
        thePeaks$symbol[is.na(thePeaks$symbol)] = ""
    }
  #  inIntron <- entrez_id
  #  inIntron[thePeaks$entrez_id != "" & thePeaks$inExon == ""] = TRUE
  #  inIntron[thePeaks$entrez_id == "" | thePeaks$inExon == TRUE] = ""
  #  thePeaks <- cbind(thePeaks, inIntron = inIntron)
    if (length(grep("FBgn", entrez_id[1])) > 0 && !missing(orgAnn) && 
        class(orgAnn) == "AnnDbBimap") {
        temp.id <- thePeaks$entrez_id
        thePeaks$entrez_id <- thePeaks$symbol
        thePeaks$symbol <- temp.id
        rm(temp.id)
    }
    unique(data.frame(thePeaks))
}

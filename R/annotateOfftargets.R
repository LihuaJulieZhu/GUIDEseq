#' Annotate offtargets with gene name
#'
#' Annotate offtargets with gene name and whether it is inside an exon
#'
#' %% ~~ If necessary, more details than the description above ~~
#'
#' @param thePeaks Output from offTargetAnalysisOfPeakRegions
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
#' @return A data frame and a tab-delimited file offTargetsInPeakRegions.xls,
#' containing all input offtargets with potential gRNA binding sites, mismatch
#' number and positions, alignment to the input gRNA and predicted cleavage
#' score, and whether the offtargets are inside an exon and associated gene
#' name.
#' @author Lihua Julie Zhu
#' @seealso GUIDEseqAnalysis
#' @references %% ~put references to the literature/web site here ~
#' @keywords utilities
#' @examples
#'
#' if (!interactive()) {
#'     library("BSgenome.Hsapiens.UCSC.hg19")
#'     library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'     library(org.Hs.eg.db)
#'     peaks <- system.file("extdata", "T2plus100OffTargets.bed",
#'         package = "CRISPRseek")
#'     gRNAs <- system.file("extdata", "T2.fa",
#'         package = "CRISPRseek")
#'     outputDir = getwd()
#'     offTargets <- offTargetAnalysisOfPeakRegions(gRNA = gRNAs, peaks = peaks,
#'         format=c("fasta", "bed"),
#'         peaks.withHeader = TRUE, BSgenomeName = Hsapiens,
#'         upstream = 20L, downstream = 20L, PAM.size = 3L, gRNA.size = 20L,
#'         orderOfftargetsBy = "predicted_cleavage_score",
#'         PAM = "NGG", PAM.pattern = "(NGG|NAG|NGA)$", max.mismatch = 2L,
#'         outputDir = outputDir,
#'         allowed.mismatch.PAM = 3, overwrite = TRUE)
#'     annotatedOfftargets <- annotateOffTargets(offTargets,
#'        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'        orgAnn = org.Hs.egSYMBOL)
#' }
#'
#' @importFrom BiocGenerics toTable
#' @importFrom GenomeInfoDb seqlevelsStyle seqnames seqlevels
#' "seqlevels<-" "seqinfo<-" "seqlevelsStyle<-"
#' @importFrom GenomicFeatures exons genes
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom S4Vectors queryHits Rle subjectHits
#'
#' @export annotateOffTargets
annotateOffTargets <-
function (thePeaks, txdb, orgAnn)
{
    thePeaks <- subset(thePeaks, !is.na(offTarget_Start) & offTarget_Start != "")
    peaks.RD <- GRanges(seqnames = Rle(thePeaks$chromosome),
        ranges = IRanges(start = as.numeric(thePeaks$offTarget_Start),
        end = as.numeric(thePeaks$offTarget_End), names = thePeaks$offTarget))
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

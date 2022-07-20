#' Offtarget Analysis of GUIDE-seq peaks
#'
#' Finding offtargets around peaks from GUIDE-seq or around any given genomic
#' regions
#'
#'
#' @param gRNA gRNA input file path or a DNAStringSet object that contains gRNA
#' plus PAM sequences used for genome editing
#' @param peaks peak input file path or a GenomicRanges object that contains
#' genomic regions to be searched for potential offtargets
#' @param format Format of the gRNA and peak input file. Currently, fasta and
#' bed are supported for gRNA and peak input file respectively
#' @param peaks.withHeader Indicate whether the peak input file contains
#' header, default FALSE
#' @param PAM.size PAM length, default 3
#' @param gRNA.size The size of the gRNA, default 20
#' @param PAM PAM sequence after the gRNA, default NGG
#' @param BSgenomeName BSgenome object. Please refer to available.genomes in
#' BSgenome package. For example, BSgenome.Hsapiens.UCSC.hg19 for hg19,
#' BSgenome.Mmusculus.UCSC.mm10 for mm10, BSgenome.Celegans.UCSC.ce6 for ce6,
#' BSgenome.Rnorvegicus.UCSC.rn5 for rn5, BSgenome.Drerio.UCSC.danRer7 for Zv9,
#' and BSgenome.Dmelanogaster.UCSC.dm3 for dm3
#' @param overlap.gRNA.positions The required overlap positions of gRNA and
#' restriction enzyme cut site, default 17 and 18 for SpCas9.
#' @param max.mismatch Maximum mismatch allowed in off target search, default 6
#' @param PAM.pattern Regular expression of protospacer-adjacent motif (PAM),
#' default to any NNN$. Set it to (NAG|NGG|NGA)$ if only outputs offtargets
#' with NAG, NGA or NGG PAM
#' @param allowed.mismatch.PAM Number of degenerative bases in the PAM.pattern
#' sequence, default to 2
#' @param outputDir the directory where the off target analysis and reports
#' will be written to
#' @param upstream upstream offset from the peak start to search for off
#' targets, default 20
#' @param downstream downstream offset from the peak end to search for off
#' targets, default 20
#' @param overwrite overwrite the existing files in the output directory or
#' not, default FALSE
#' @param weights a numeric vector size of gRNA length, default c(0, 0, 0.014,
#' 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732,
#' 0.828, 0.615, 0.804, 0.685, 0.583) for SPcas9 system, which is used in Hsu
#' et al., 2013 cited in the reference section. Please make sure that the
#' number of elements in this vector is the same as the gRNA.size, e.g., pad 0s
#' at the beginning of the vector.
#' @param orderOfftargetsBy criteria to order the offtargets by and the top one
#' will be kept if keepTopOfftargetsOnly is set to TRUE.  If set to
#' predicted_cleavage_score (descending order), the offtarget with the highest
#' predicted cleavage score for each peak will be kept.  If set to n.mismatch
#' (ascending order), the offtarget with the smallest number of mismatch to the
#' target sequence for each peak will be kept.
#' @param descending No longer used. In the descending or ascending order.
#' Default to order by predicted cleavage score in descending order and number
#' of mismatch in ascending order When altering orderOfftargetsBy order, please
#' also modify descending accordingly
#' @param keepTopOfftargetsOnly Output all offtargets or the top offtarget per
#' peak using the orderOfftargetsBy criteria, default to the top offtarget
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
#' possible types of single nucleotide mismatch at each position of the gRNA.
#' By default, using the supplemental Table 19 from Doench et al., Nature
#' Biotechnology 2016
#' @param n.cores.max Indicating maximum number of cores to use in multi core
#' mode, i.e., parallel processing, default 1 to disable multicore processing
#' for small dataset.
#' @return a tab-delimited file offTargetsInPeakRegions.tsv, containing all
#' input peaks with potential gRNA binding sites, mismatch number and
#' positions, alignment to the input gRNA and predicted cleavage score.
#' @author Lihua Julie Zhu
#' @seealso GUIDEseq
#' @references Patrick D Hsu, David A Scott, Joshua A Weinstein, F Ann Ran,
#' Silvana Konermann, Vineeta Agarwala, Yinqing Li, Eli J Fine, Xuebing Wu,
#' Ophir Shalem,Thomas J Cradick, Luciano A Marraffini, Gang Bao & Feng Zhang
#' (2013) DNA targeting specificity of rNA-guided Cas9 nucleases.  Nature
#' Biotechnology 31:827-834 Lihua Julie Zhu, Benjamin R. Holmes, Neil Aronin
#' and Michael Brodsky.  CRISPRseek: a Bioconductor package to identify
#' target-specific guide RNAs for CRISPR-Cas9 genome-editing systems. Plos One
#' Sept 23rd 2014 Lihua Julie Zhu (2015). Overview of guide RNA design tools
#' for CRISPR-Cas9 genome editing technology. Frontiers in Biology August 2015,
#' Volume 10, Issue 4, pp 289-296
#' @keywords misc
#' @examples
#'
#' #### the following example is also part of annotateOffTargets.Rd
#' if (interactive())
#' {
#'     library("BSgenome.Hsapiens.UCSC.hg19")
#'     library(GUIDEseq)
#'     peaks <- system.file("extdata", "T2plus100OffTargets.bed",
#'         package = "CRISPRseek")
#'     gRNAs <- system.file("extdata", "T2.fa",
#'         package = "CRISPRseek")
#'     outputDir = getwd()
#'     offTargets <- offTargetAnalysisOfPeakRegions(gRNA = gRNAs, peaks = peaks,
#'         format=c("fasta", "bed"),
#'         peaks.withHeader = TRUE, BSgenomeName = Hsapiens,
#'         upstream = 25L, downstream = 25L, PAM.size = 3L, gRNA.size = 20L,
#'         orderOfftargetsBy = "predicted_cleavage_score",
#'         PAM = "NGG", PAM.pattern = "(NGG|NAG|NGA)$", max.mismatch = 2L,
#'         outputDir = outputDir,
#'         allowed.mismatch.PAM = 3, overwrite = TRUE
#'    )
#' }
#' @importFrom utils read.table write.table
#' @importFrom IRanges IRanges nearest
#' @importFrom hash hash
#' @importFrom CRISPRseek compare2Sequences
#' @importFrom parallel detectCores makeCluster parLapply
#' stopCluster
#' @export offTargetAnalysisOfPeakRegions
offTargetAnalysisOfPeakRegions <-
    function(gRNA, peaks,
    format=c("fasta", "bed"),
    peaks.withHeader = FALSE,
    BSgenomeName,
    overlap.gRNA.positions = c(17,18),
    upstream = 25L,
    downstream =25L,
    PAM.size = 3L,
    gRNA.size = 20L,
    PAM = "NGG",
    PAM.pattern = "NNN$",
    max.mismatch = 6L,
    outputDir,
    allowed.mismatch.PAM = 2L,
    overwrite = TRUE,
    weights = c(0, 0, 0.014, 0, 0, 0.395,
    0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,
    0.804, 0.685, 0.583),
    orderOfftargetsBy = c("predicted_cleavage_score", "n.mismatch"),
    descending = TRUE,
    keepTopOfftargetsOnly = TRUE,
    scoring.method = c("Hsu-Zhang", "CFDscore"),
        subPAM.activity = hash( AA = 0, AC = 0, AG = 0.259259259, AT = 0,
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
     n.cores.max = 1
   )
{
    orderOfftargetsBy <- match.arg(orderOfftargetsBy)
    thePeaks <- read.table(peaks, sep="\t", header = peaks.withHeader,
        stringsAsFactors = FALSE)
    if (format[2] == "bed")
    {
       if (dim(thePeaks)[2] == 3)
       {
	  pnames <- paste(thePeaks[,1], thePeaks[,2], sep = ":")
	  thePeaks <- cbind(thePeaks, pnames, pnames, pnames)
          thePeaks[, 5] <-  0
          thePeaks[, 6] <- "+"
       }
       if (dim(thePeaks)[2] >= 4)
       {
           colnames(thePeaks)[1:4] <- c("chromosome", "peak_start", "peak_end",
                "names")
           if (dim(thePeaks)[2] == 4)
           {
		thePeaks <- cbind(thePeaks, thePeaks[, 4], thePeak[, 5])
                thePeaks[, 5] <-  0
                thePeaks[, 6] <- "+"
	   }
       }
        if(dim(thePeaks)[2] >= 5)
        {
           colnames(thePeaks)[5] = "peak_score"
           if (dim(thePeaks)[2] == 5)
           {
                thePeaks <- cbind(thePeaks, thePeaks[, 5])
                thePeaks[, 6] <- "+"
           }
        }
        if(dim(thePeaks)[2] >= 6)
            colnames(thePeaks)[6] = "peak_strand"
    }
    else
    {
        stop("only bed file with at least 4 columns are supported")
    }
    thePeaks[, 4] <- gsub(" ", "", thePeaks[,4])
    TS2 <- tryCatch(
        (compare2Sequences(inputFile1Path = gRNA, inputFile2Path = peaks,
        findgRNAsWithREcutOnly = FALSE, format = format,
        header = peaks.withHeader,
        BSgenomeName = BSgenomeName,
        upstream = upstream,
        downstream = downstream,
        minREpatternSize = 4,
        findgRNAs = c(FALSE, FALSE),
        removegRNADetails = c(TRUE, FALSE),
        exportAllgRNAs =  "no",
        annotatePaired = FALSE,
        searchDirection = "1to2",
        overlap.gRNA.positions = overlap.gRNA.positions,
        findPairedgRNAOnly = FALSE,
        min.gap = 0, max.gap = 20, gRNA.name.prefix = "gRNA",
        PAM.size = PAM.size,
        gRNA.size = gRNA.size, PAM = PAM, PAM.pattern = PAM.pattern,
        max.mismatch = max.mismatch,
        outputDir = outputDir, foldgRNAs = FALSE,
        allowed.mismatch.PAM = allowed.mismatch.PAM, overwrite = overwrite,
        weights = weights,
        scoring.method = scoring.method,
        subPAM.activity = subPAM.activity,
        subPAM.position = subPAM.position,
        PAM.location = PAM.location,
        mismatch.activity.file = mismatch.activity.file
        )),
        error = function(e) {
            message(e)
            return(NA)
        })

    if ( exists("TS2") && length(TS2) > 1)
    {
        if (dim(TS2)[1] > 200)
            n.cores <- min(n.cores.max, floor(detectCores() /3) + 1)
        else
            n.cores <- 1
        if (n.cores > 1)
        {
            cl <- makeCluster(n.cores)

            names <- as.character(unlist(parLapply(cl, as.character(TS2$offTarget),
                function(temp) {
                    temp1 <- strsplit(temp,":")[[1]]
                    end.ind <- length(temp1) - 1
                    paste(temp1[1:end.ind], collapse=":")
            })))
           TS2 <- cbind(names = names, TS2)
           excluding.columns = which(colnames(TS2) %in%
               c("scoreForSeq1", "targetInSeq1", "gRNAefficacy", "scoreDiff"))
           TS2 <- TS2[, -excluding.columns]
           colnames(TS2)[colnames(TS2) == "scoreForSeq2"] = "predicted_cleavage_score"
           colnames(TS2)[colnames(TS2) == "targetInSeq2"] = "offTarget_sequence"
           offTargetOffset <- do.call(rbind, parLapply(cl, as.character(TS2$offTarget),
               function(temp)
               {
               temp1 <- strsplit(temp,":")[[1]]
               start.ind <- length(temp1)
               as.numeric(strsplit(temp1[start.ind], "-")[[1]][1:2])
             }))
            stopCluster(cl)
        }
        else
        {
            names <- as.character(unlist(lapply(as.character(TS2$offTarget),
                function(temp) {
                    temp1 <- strsplit(temp,":")[[1]]
                    end.ind <- length(temp1) - 1
                    paste(temp1[1:end.ind], collapse=":")
            })))
            TS2 <- cbind(names = names, TS2)
            excluding.columns = which(colnames(TS2) %in%
               c("scoreForSeq1", "targetInSeq1", "gRNAefficacy", "scoreDiff"))
            TS2 <- TS2[, -excluding.columns]
            colnames(TS2)[colnames(TS2) == "scoreForSeq2"] = "predicted_cleavage_score"
            colnames(TS2)[colnames(TS2) == "targetInSeq2"] = "offTarget_sequence"
            offTargetOffset <- do.call(rbind, lapply(as.character(TS2$offTarget),
               function(temp)
               {
               temp1 <- strsplit(temp,":")[[1]]
               start.ind <- length(temp1)
               as.numeric(strsplit(temp1[start.ind], "-")[[1]][1:2])
             }))
        }
        TS2 <- cbind(TS2, offTarget_Start = offTargetOffset[,1],
           offTarget_End = offTargetOffset[,2])
        offtargets <- merge(TS2, thePeaks, by = "names", all = TRUE)
        offtargets$offTarget_Start <-
            ifelse(as.character(offtargets$peak_strand) == "+",
            offtargets$offTarget_Start - upstream + offtargets$peak_start - 1,
            upstream - offtargets$offTarget_End + offtargets$peak_end + 1)
        offtargets$offTarget_End <-
            offtargets$offTarget_Start + PAM.size  + gRNA.size  - 1

        offtargets.minus.minus <-
            subset(offtargets, as.character(offtargets$peak_strand) == "-" &
            as.character(offtargets$offTargetStrand) == "-")
        offtargets.minus.plus <- subset(offtargets,
            as.character(offtargets$peak_strand) == "-" &
            as.character(offtargets$offTargetStrand) == "+")
        if (dim(offtargets.minus.minus)[1] > 0)
            offtargets.minus.minus$offTargetStrand <- "+"
        if (dim(offtargets.minus.plus)[1] > 0)
            offtargets.minus.plus$offTargetStrand<- "-"

        offtargets <- rbind(subset(offtargets,
            as.character(offtargets$peak_strand) == "+" |
            as.character(offtargets$peak_strand) == "*" |
            is.na(offtargets$offTargetStrand)),
            offtargets.minus.minus, offtargets.minus.plus)
        if (keepTopOfftargetsOnly)
        {
            peaks.without.offtargets <- subset(offtargets, is.na(gRNAPlusPAM) | gRNAPlusPAM == "")
            offtargets <- subset(offtargets, !is.na(gRNAPlusPAM) & gRNAPlusPAM != "")
            if (dim(offtargets)[1] > 1)
            {
                offtargets$predicted_cleavage_score <-
                    as.numeric(as.character(offtargets$predicted_cleavage_score))
                offtargets$n.mismatch <- as.numeric(as.character(offtargets$n.mismatch))
                temp <- aggregate(predicted_cleavage_score ~ names, offtargets, max)
                if (orderOfftargetsBy == "n.mismatch")
                     temp <- aggregate(n.mismatch ~ names, offtargets, min)
                offtargets <- merge(offtargets, temp)
             }
            ###### keep only one nearest offtarget for each peak
             temp <- as.data.frame(table(offtargets$names))
             names.notUnique <- subset(offtargets,
                  offtargets$names %in% temp[temp[,2] > 1, 1])
             offtargets <- subset(offtargets,
                  offtargets$names %in% temp[temp[,2] ==1, 1])
             if (dim(temp[temp[,2] > 1, ])[1] > 0)
             {
                 for (nu.name in temp[temp[,2] >1, 1])
                 {
                      notUnique <- subset(names.notUnique, names == nu.name)
                      this.peak <- IRanges(start = notUnique$peak_start[1],
                          end = notUnique$peak_end[1])
                      nu.ots <- IRanges(start = notUnique$offTarget_Start,
                          end = notUnique$offTarget_End)
                      offtargets <- rbind( offtargets,
                          notUnique[nearest(this.peak, nu.ots),])
                 }
             }
           ###### keep only one nearest peak for each offtarget
             offtargets$offTarget <-
                 paste(as.character(offtargets$chromosome),
                     as.character(offtargets$offTargetStrand),
                     as.character(offtargets$offTarget_Start),
                     as.character(offtargets$offTarget_End), sep=":")
             temp <- aggregate(peak_score ~ offTarget, offtargets, max)
             offtargets <- merge(offtargets, temp)
             if (dim(offtargets)[1] > dim(temp)[1])
             {
	         temp <- as.data.frame(table(offtargets$offTarget))
                 offtargets.notUnique <- subset(offtargets,
                     offtargets$offTarget %in% temp[temp[,2] > 1, 1])
                 offtargets <- subset(offtargets,
                     offtargets$offTarget %in% temp[temp[,2] ==1, 1])
                 for (ot in temp[temp[,2] >1, 1])
                 {
                     notUnique <- subset(offtargets.notUnique, offTarget == ot)
                     this.ot <- IRanges(start = notUnique$offTarget_Start[1],
                          end = notUnique$offTarget_End[1])
                     nu.peaks <- IRanges(start = notUnique$peak_start, end = notUnique$peak_end)
		     offtargets <- rbind( offtargets,
                        notUnique[nearest(this.ot, nu.peaks),])
                 }
             }

             #offtargets <- rbind(offtargets, peaks.without.offtargets)
        }
        write.table(offtargets,
            file = file.path(outputDir,"offTargetsInPeakRegions.xls"),
            sep="\t", row.names = FALSE)
        unlink(file.path(outputDir, "scoresFor2InputSequences.xls"))
        offtargets
     }
    else
    {
        write.table(thePeaks,
            file = file.path(outputDir,"offTargetsInPeakRegions.xls"),
            sep="\t", row.names = FALSE)
        thePeaks
    }
}

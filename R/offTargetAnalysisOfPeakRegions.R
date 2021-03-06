offTargetAnalysisOfPeakRegions <-
    function(gRNA, peaks,
    format=c("fasta", "bed"),
    peaks.withHeader = FALSE,
    BSgenomeName,
    overlap.gRNA.positions = c(17,18),
    upstream = 20L,
    downstream =20L,
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

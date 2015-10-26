getPeaks <-
function(gr, window.size = 20L, step = 20L, bg.window.size = 5000L,
    min.reads = 10L, min.SNratio = 2, maxP = 0.05,
    n.cores.max = 6, stats = c("poisson", "nbinom"),
    p.adjust.methods = 
    c( "none", "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr"))
{
    message("Validating input ...\n");
    if (missing(gr)) {
        stop("Missing required argument gr!")
    }
    if (class(gr) != "GRanges" ) {
        stop("No valid gr passed in. It needs to be GRanges object")
    }
    stats <- match.arg(stats)
    p.adjust.methods <- match.arg(p.adjust.methods)
    plus.gr = subset(gr, strand(gr) == "+")
    minus.gr = subset(gr, strand(gr) == "-")
    if (length(plus.gr) >= 2)
    {
        message("parpare for plus strand ...\n");
        plus.runsum <- .getStrandedCoverage(plus.gr, window.size = window.size,
            bg.window.size = bg.window.size, min.reads = min.reads,
            strand = "+")
        if (length(minus.gr) >= 2)
        {
            message("prepare for minus strand ...\n");
            minus.runsum <- .getStrandedCoverage(minus.gr,
                window.size = window.size,
                bg.window.size = bg.window.size, min.reads = min.reads,
                strand = "-")
            both.runsum <- c(plus.runsum, minus.runsum)
        }
        else
        {
            both.runsum <- plus.runsum
        }
    }
    else if (length(minus.gr) >= 2)
    {
        message("prepare for minus strand ...\n");
        both.runsum <- .getStrandedCoverage(minus.gr, window.size = window.size,
            bg.window.size = bg.window.size, 
            min.reads = min.reads, strand = "-")
    }
    else
        stop("too few reads for peak calling! \n")
    message("call peaks ...\n");
    if (!exists("both.runsum"))
        stop("no peaks found! \n")
    if (stats == "poisson")
    {
        both.runsum$p.value <- ppois(as.numeric(both.runsum$count),
            lambda = as.numeric(both.runsum$bg), lower.tail = FALSE,
            log.p = FALSE)
    }
    else if (stats == "nbinom")
    {
        both.runsum$p.value <- pnbinom(as.numeric(both.runsum$count),
            mu = as.numeric(both.runsum$bg),
            size = window.size, lower.tail = FALSE, log.p = FALSE)
    }
    both.runsum$SNratio <- 
        as.numeric(both.runsum$count)/as.numeric(both.runsum$bg)
    n.cores <- detectCores()
    n.cores <- min(n.cores.max, n.cores, 
       length(seqnames(seqinfo(both.runsum))))
    if (n.cores > 1)
    {
        cl <- makeCluster(n.cores)
        clusterExport(cl, varlist = c("both.runsum", ".locMaxPos"),
            envir = environment()) 
        clusterExport(cl, varlist =  c("window.size", "step", "min.reads"),
            envir = environment())
        local.max.gr <- do.call(c, parLapply(cl, seqnames(seqinfo(both.runsum)),
            function(chr) {
            message("processing chromosome", chr, "\n")
            this.gr <- subset(both.runsum, seqnames(both.runsum) == chr & 
                strand(both.runsum) == "+")
            minus.gr <- subset(both.runsum, seqnames(both.runsum) == chr & 
                strand(both.runsum) == "-")
            #max.pos <- which(diff(sign(diff(as.numeric(
                #as.character(this.gr$count)))))==-2)+1
            if (length(this.gr) >= 1) {
                max.pos <- .locMaxPos(this.gr, window.size = window.size,
                    step = step, min.reads = min.reads)
                plus.grs <- this.gr[start(this.gr) %in% max.pos]
                max.pos2 <- .locMaxPos(plus.grs, window.size = window.size,
                    step = step, min.reads = min.reads)
                plus.grs2 <- plus.grs[start(plus.grs) %in% max.pos2]
                if (length(minus.gr) >= 1) {
                    max.pos <- .locMaxPos(minus.gr, window.size = window.size,
                        step = step, min.reads = min.reads)
                    minus.grs <- minus.gr[start(minus.gr) %in% max.pos]
                    max.pos2 <- .locMaxPos(minus.grs, window.size = window.size,
                        step = step, min.reads = min.reads)
                    minus.grs2 <- minus.grs[start(minus.grs) %in% max.pos2]
                    max.gr <- c(plus.grs2, minus.grs2)
                }
                else {
                    max.gr <- plus.grs2
                }
            }
            else if (length(minus.gr) >= 1) {
                max.pos <- .locMaxPos(minus.gr, window.size = window.size,
                    step = step, min.reads = min.reads)
                minus.grs <- minus.gr[start(minus.gr) %in% max.pos]
                max.pos2 <- .locMaxPos(minus.grs, window.size = window.size,
                    step = step, min.reads = min.reads)
                max.gr <- minus.grs[start(minus.grs) %in% max.pos2]
            }
            else {
                max.gr = GRanges()
            }
            max.gr
        }))
        stopCluster(cl)
    }
    else
    {
        local.max.gr <- do.call(c, lapply(seqnames(seqinfo(both.runsum)),
            function(chr) {
            message("processing chromosome", chr, "\n")
            this.gr <- both.runsum[seqnames(both.runsum) == chr &
                strand(both.runsum) == "+"]
            minus.gr <- both.runsum[seqnames(both.runsum) == chr &
                strand(both.runsum) == "-"]
            #max.pos <- which(diff(sign(diff(as.numeric(
                #as.character(this.gr$count)))))==-2)+1
            if (length(this.gr) >= 1) {
                max.pos <- .locMaxPos(this.gr, window.size = window.size,
                    step = step, min.reads = min.reads)
                plus.grs2 <- this.gr[start(this.gr) %in% max.pos] 
                #plus.grs <- this.gr[start(this.gr) %in% max.pos]
                #max.pos2 <- .locMaxPos(plus.grs, window.size = window.size,
                    #step = step, min.reads = min.reads)
                #plus.grs2 <- plus.grs[start(plus.grs) %in% max.pos2]
                if (length(minus.gr) >= 1) {
                    max.pos <- .locMaxPos(minus.gr, window.size = window.size,
                        step = step, min.reads = min.reads)
                    minus.grs2 <- minus.gr[start(minus.gr) %in% max.pos]
                    #minus.grs <- minus.gr[start(minus.gr) %in% max.pos]
                    #max.pos2 <- .locMaxPos(minus.grs, window.size = window.size,
                    #    step = step, min.reads = min.reads)
                    #minus.grs2 <- minus.grs[start(minus.grs) %in% max.pos2]
                    max.gr <- c(plus.grs2, minus.grs2)
                }
                else {
                    max.gr <- plus.grs2
                }
            }
            else if (length(minus.gr) >= 1) {
                max.pos <- .locMaxPos(minus.gr, window.size = window.size,
                    step = step, min.reads = min.reads)
                minus.grs <- minus.gr[start(minus.gr) %in% max.pos]
                #max.pos2 <- .locMaxPos(minus.grs, window.size = window.size,
                #    step = step, min.reads = min.reads)
                #max.gr <- minus.grs[start(minus.grs) %in% max.pos2]
                max.gr <- minus.grs[start(minus.grs) %in% max.pos]
            }
            else {
                max.gr = GRanges()
            }
            max.gr
        }))
    }
    both.runsum.bk <- both.runsum
    both.runsum <- local.max.gr
   
    if (length(both.runsum) > 0)
    {
        if (p.adjust.methods != "none")
        {
            both.runsum$adjusted.p.value <- p.adjust(both.runsum$p.value,
                method = p.adjust.methods)
            peaks <- subset(both.runsum, both.runsum$adjusted.p.value <= maxP & 
                both.runsum$SNratio >= min.SNratio)
        }
        else
        {
            peaks <- subset(both.runsum, both.runsum$p.value <= maxP & 
                both.runsum$SNratio >= min.SNratio)
        }
    }
    else {
       peaks <- both.runsum
    }
    list(peaks = peaks, both.runsum.bk = both.runsum.bk, 
        summarized.count = both.runsum)
}

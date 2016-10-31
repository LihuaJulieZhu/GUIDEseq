createBarcodeFasta <- function(p5.index, p7.index,
    reverse.p7 = TRUE, reverse.p5 = FALSE,
    header = FALSE, outputFile = "barcodes.fa")
{
	p7 <- read.table(p7.index, header = header, stringsAsFactors=FALSE)
	if (reverse.p7)
           {
	    p7.rev <- unlist(lapply(p7[,1], function(x) 
                  {as.character(reverseComplement(DNAString(x)))}
	    ))
            }
            else
            {
	     p7.rev <- p7[,1]
             }
	p5 <- read.table(p5.index, header = header, stringsAsFactors=FALSE)
           if (reverse.p5)
           {
	    p5.rev <- unlist(lapply(p5[,1], function(x) 
                  {as.character(reverseComplement(DNAString(x)))}
	    ))
            }
            else
            {
	     p5.rev <- p5[,1]
             }
	n.barcode <- length(p5.rev) * length(p7.rev)
	barcode <- character(n.barcode)
	k <- 0;
        for (i in p7.rev) 
        {
            for (j in p5.rev) 
            {
		k <- k +1
                barcode[k] <-paste(i, j, sep="")
             }
         }
        seqs <- cbind(barcode, sequence = as.character(barcode))
        seqs[,1] = paste(">", seqs[,1], sep="")
         write.table(matrix(t(seqs), ncol =1, nrow = dim(seqs)[1] * 2), 
             file = outputFile, sep="\n", row.names=FALSE, 
             col.names=FALSE, quote=FALSE)
}

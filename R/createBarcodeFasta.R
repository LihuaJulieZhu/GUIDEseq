#' Create barcode as fasta file format for building bowtie1 index
#' 
#' Create barcode as fasta file format for building bowtie1 index to assign
#' reads to each library with different barcodes. The bowtie1 index has been
#' built for the standard GUIDE-seq protocol using the standard p5 and p7
#' index. It can be downloaded at
#' http://mccb.umassmed.edu/GUIDE-seq/barcode.bowtie1.index.tar.gz
#' 
#' 
#' @param p5.index A text file with one p5 index sequence per line
#' @param p7.index A text file with one p7 index sequence per line
#' @param header Indicate whether there is a header in the p5.index and
#' p7.index files.  Default to FALSE
#' @param reverse.p7 Indicate whether to reverse p7 index, default to TRUE for
#' standard GUIDE-seq experiments
#' @param reverse.p5 Indicate whether to reverse p5 index, default to FALSE for
#' standard GUIDE-seq experiments
#' @param outputFile Give a name to the output file where the generated
#' barcodes are written. This file can be used to build bowtie1 index for
#' binning reads.
#' @note Create barcode file to be used to bin the reads sequenced in a pooled
#' lane
#' @author Lihua Julie Zhu
#' @references %% ~put references to the literature/web site here ~
#' @keywords manip utilities
#' @examples
#' 
#'     p7 <- system.file("extdata", "p7.index",
#'            package = "GUIDEseq") 
#'     p5 <- system.file("extdata", "p5.index",
#'            package = "GUIDEseq")
#'     outputFile <- "barcodes.fa" 
#'     createBarcodeFasta(p5.index = p5, p7.index = p7, reverse.p7 = TRUE,
#'         reverse.p5 = FALSE, outputFile = outputFile)
#' 
#' @export createBarcodeFasta
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

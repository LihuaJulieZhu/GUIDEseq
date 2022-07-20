#' Create barcodes from the p5 and p7 index used for each sequencing lane
#'
#' Create barcodes from the p5 and p7 index for assigning reads to each barcode
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
#' barcodes are written
#' @return DNAStringSet
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
#'     outputFile <- "usedBarcode"
#'     getUsedBarcodes(p5.index = p5, p7.index = p7, reverse.p7 = TRUE,
#'         reverse.p5 = FALSE, outputFile = outputFile)
#' @importFrom Biostrings DNAString DNAStringSet reverseComplement
#' @importFrom utils read.table write.table
#' @export getUsedBarcodes
#'
getUsedBarcodes <- function(p5.index, p7.index, header = FALSE,
   reverse.p7 = TRUE, reverse.p5 = FALSE, outputFile)
{
   if (file.exists(p5.index))
   {
      p5 <- read.table(p5.index, header = header, stringsAsFactors=FALSE)
      p5 <- p5[,1]
   }
   else
   {
      p5 <- p5.index
   }
   if (file.exists(p7.index))
   {
      p7 <- read.table(p7.index, header = header, stringsAsFactors=FALSE)
      p7 <- p7[,1]
   }
   else
   {
      p7 <- p7.index
   }
   if (reverse.p7)
   {
      p7.rev <- unlist(lapply(p7, function(x)
         {as.character(reverseComplement(DNAString(x)))}
         ))
    }
    else
    {
       p7.rev <- p7
    }
    if (reverse.p5)
    {
       p5.rev <- unlist(lapply(p5, function(x)
          {as.character(reverseComplement(DNAString(x)))}
	  ))
    }
   else
   {
      p5.rev <- p5
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
   write.table(paste(barcode, collapse = " "), quote = FALSE,
     file = outputFile, row.names = FALSE, col.names = FALSE)
   DNAStringSet(barcode)
}

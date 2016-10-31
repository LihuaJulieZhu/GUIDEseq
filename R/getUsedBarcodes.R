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

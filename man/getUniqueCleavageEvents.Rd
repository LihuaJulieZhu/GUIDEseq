\name{getUniqueCleavageEvents}
\alias{getUniqueCleavageEvents}
\title{
Using UMI sequence to obtain the starting sequence library
}
\description{
PCR amplification often leads to biased representation of the starting sequence
population. To track the sequence tags present in the initial sequence library, 
a unique molecular identifier (UMI) is added to the 5 prime of each sequence in 
the staring library. This function uses the UMI sequence plus the first few 
sequence from R1 reads to obtain the starting sequence library.
}
\usage{
getUniqueCleavageEvents(alignment.inputfile, umi.inputfile, 
 alignment.format = c("auto", "bam", "bed"),
 umi.header = FALSE, read.ID.col = 1, 
umi.col = 2, umi.sep = "\t", keep.chrM = FALSE,
keep.R1only = TRUE, keep.R2only = TRUE, 
concordant.strand = TRUE, max.paired.distance = 1000, 
min.mapping.quality = 30, max.R1.len = 130, max.R2.len = 130, 
apply.both.max.len = FALSE, same.chromosome = TRUE, 
distance.inter.chrom = -1, min.R1.mapped = 20, min.R2.mapped = 20, 
apply.both.min.mapped = FALSE, max.duplicate.distance = 0, 
umi.plus.R1start.unique = TRUE, umi.plus.R2start.unique = TRUE,
n.cores.max = 6)
}
\arguments{
  \item{alignment.inputfile}{
The alignment file. Currently supports bed output file with CIGAR information.
Suggest run the workflow binReads.sh, which sequentially runs barcode binning,
adaptor removal, alignment to genome, alignment quality filtering, and bed file
conversion. Please download the workflow function and its helper scripts at
http://mccb.umassmed.edu/GUIDE-seq/binReads/
}
  \item{umi.inputfile}{
A text file containing at least two columns, one is the read identifier and the
other is the UMI or UMI plus the first few bases of R1 reads. Suggest use 
getUMI.sh to generate this file. Please download the script and its helper 
scripts at http://mccb.umassmed.edu/GUIDE-seq/getUMI/
}
  \item{alignment.format}{
The format of the alignment input file. Currently only
bam and bed file format is supported. BED format will be deprecated soon.
}
  \item{umi.header}{
Indicates whether the umi input file contains a header line or not. 
Default to FALSE
}
  \item{read.ID.col}{
The index of the column containing the read identifier in the umi input file, 
default to 1
}
  \item{umi.col}{
The index of the column containing the umi or umi plus the first few bases 
of sequence from the R1 reads, default to 2
}
  \item{umi.sep}{
column separator in the umi input file, default to tab
}
\item{keep.chrM}{
Specify whether to include alignment from chrM.
Default FALSE 
}
  \item{keep.R1only}{
Specify whether to include alignment with only R1 without paired R2. 
Default TRUE
}
  \item{keep.R2only}{
Specify whether to include alignment with only R2 without paired R1. 
Default TRUE
}
  \item{concordant.strand}{
Specify whether the R1 and R2 should be aligned to the same strand or opposite
strand. Default opposite.strand (TRUE) 
}
  \item{max.paired.distance}{
Specify the maximum distance allowed between paired R1 and R2 reads. 
Default 1000 bp
}
  \item{min.mapping.quality}{
Specify min.mapping.quality of acceptable alignments
}
  \item{max.R1.len}{
The maximum retained R1 length to be considered for downstream analysis, 
default 130 bp. Please note that default of 130 works well when the read 
length 150 bp. Please also note that retained R1 length is not necessarily 
equal to the mapped R1 length
}
  \item{max.R2.len}{
The maximum retained R2 length to be considered for downstream analysis, 
default 130 bp. Please note that default of 130 works well when the read 
length 150 bp. Please also note that retained R2 length is not necessarily 
equal to the mapped R2 length
}
  \item{apply.both.max.len}{
Specify whether to apply maximum length requirement to both R1 and R2 reads,
default FALSE
}
  \item{same.chromosome}{
Specify whether the paired reads are required to align to the same chromosome,
default TRUE
}
  \item{distance.inter.chrom}{
Specify the distance value to assign to the paired reads that are aligned to
different chromosome, default -1
}
  \item{min.R1.mapped}{
The maximum mapped R1 length to be considered for downstream analysis, 
default 30 bp. 
}
  \item{min.R2.mapped}{
The maximum mapped R2 length to be considered for downstream analysis, 
default 30 bp.
}
  \item{apply.both.min.mapped}{
Specify whether to apply minimum mapped length requirement to both R1 and 
R2 reads, default FALSE
}
  \item{max.duplicate.distance}{
Specify the maximum distance apart for two reads to be considered as 
duplicates, default 0. Currently only 0 is supported
}
  \item{umi.plus.R1start.unique}{
To specify whether two mapped reads are considered as unique if both 
containing the same UMI and same alignment start for R1 read, default TRUE.
}
  \item{umi.plus.R2start.unique}{
To specify whether two mapped reads are considered as unique if both 
containing the same UMI and same alignment start for R2 read, default TRUE.
}
\item{n.cores.max}{Indicating maximum number of cores to use in multi core mode, 
    i.e., parallel processing, default 6. Please set it to 1 to disable multicore
    processing for small dataset.
}
}
\value{
\item{cleavage.gr }{Cleavage sites with one site per UMI as GRanges with 
metadata column total set to 1 for each range}
\item{unique.umi.plus.R2}{a data frame containing unique cleavage site 
from R2 reads mapped to plus strand with
   the following columns
   chr.y (chromosome of readSide.y/R2 read)   
   chr.x (chromosome of readSide.x/R1 read)
   strand.y (strand of readSide.y/R2 read)
   strand.x (strand of readSide.x/R1 read)
   start.y (start of readSide.y/R2 read)
   end.x (start of readSide.x/R1 read)
   UMI (unique molecular identifier (umi) or umi with the first few bases 
   of R1 read)
}
\item{unique.umi.minus.R2}{a data frame containing unique cleavage site 
from R2 reads mapped to minus strand
   with the following columns
   chr.y (chromosome of readSide.y/R2 read)   
   chr.x (chromosome of readSide.x/R1 read)
   strand.y (strand of readSide.y/R2 read)
   strand.x (strand of readSide.x/R1 read)
   end.y (end of readSide.y/R2 read)
   start.x (start of readSide.x/R1 read)
   UMI (unique molecular identifier (umi) or umi with the first few bases 
   of R1 read)
}
\item{unique.umi.plus.R1}{a data frame containing unique cleavage site from 
R1 reads mapped to minus strand
   without corresponding R2 reads mapped to the plus strand, with the 
   following columns
   chr.y (chromosome of readSide.y/R2 read)   
   chr.x (chromosome of readSide.x/R1 read)
   strand.y (strand of readSide.y/R2 read)
   strand.x (strand of readSide.x/R1 read)
   start.x (start of readSide.x/R1 read)
   start.y (start of readSide.y/R2 read)
   UMI (unique molecular identifier (umi) or umi with the first few bases 
   of R1 read)
}
\item{unique.umi.minus.R1}{a data frame containing unique cleavage site from 
R1 reads mapped to plus strand
   without corresponding R2 reads mapped to the minus strand, with the 
   following columns
   chr.y (chromosome of readSide.y/R2 read)   
   chr.x (chromosome of readSide.x/R1 read)
   strand.y (strand of readSide.y/R2 read)
   strand.x (strand of readSide.x/R1 read)
   end.x (end of readSide.x/R1 read)
   end.y (end of readSide.y/R2 read)
   UMI (unique molecular identifier (umi) or umi with the first few bases 
   of R1 read)
}
\item{all.umi}{a data frame containing all the mapped reads with the following 
columns. 
   readName (read ID),   chr.x (chromosome of readSide.x/R1 read), 
   start.x (start of eadSide.x/R1 read),
   end.x (end of eadSide.x/R1 read), 
   mapping.qual.x (mapping quality of readSide.x/R1 read),  
   strand.x (strand of readSide.x/R1 read),
   cigar.x (CIGAR of readSide.x/R1 read) , readSide.x (1/R1), 
   chr.y (chromosome of readSide.y/R2 read)        
   start.y (start of readSide.y/R2 read), 
   end.y (end of readSide.y/R2 read), 
   mapping.qual.y (mapping quality of readSide.y/R2 read),
   strand.y (strand of readSide.y/R2 read), 
   cigar.y (CIGAR of readSide.y/R2 read), readSide.y (2/R2)   
   R1.base.kept (retained R1 length), R2.base.kept (retained R2 length),
   distance (distance between mapped R1 and R2), 
   UMI (unique molecular identifier (umi) or umi with the first few bases 
   of R1 read)
}
}
\references{
Shengdar Q Tsai and J Keith Joung et al. GUIDE-seq enables genome-wide 
profiling of off-target cleavage by CRISPR-Cas nucleases. Nature Biotechnology
33, 187 to 197 (2015)
}
\author{
Lihua Julie Zhu
}
\seealso{
getPeaks
}
\examples{
    if(interactive())
    {
        umiFile <- system.file("extdata", "UMI-HEK293_site4_chr13.txt",
           package = "GUIDEseq")
        alignFile <- system.file("extdata","bowtie2.HEK293_site4_chr13.sort.bam" ,
            package = "GUIDEseq")
        cleavages <- getUniqueCleavageEvents(
            alignment.inputfile = alignFile , umi.inputfile = umiFile,
            n.cores.max = 1)
        names(cleavages)
        #output a summary of duplicate counts for sequencing saturation assessment
        table(cleavages$umi.count.summary$n)
    }
}
\keyword{ misc }

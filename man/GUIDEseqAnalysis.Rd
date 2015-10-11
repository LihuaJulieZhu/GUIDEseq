\name{GUIDEseqAnalysis}
\alias{GUIDEseqAnalysis}
\title{
Analysis pipeline for GUIDE-seq dataset
}
\description{
A wrapper function that uses the UMI sequence plus the first few bases of each
    sequence from R1 reads to estimate the starting sequence library, piles up
    reads with a user defined window and step size, identify the cleavage sites,
    merge cleavage sites from plus strand and minus strand, followed by off
    target analysis of extended regions around the identified cleavage sites.
}
\usage{
GUIDEseqAnalysis(alignment.inputfile, umi.inputfile,
    alignment.format = "bed", 
    umi.header = FALSE, read.ID.col = 1,
    umi.col = 2, umi.sep = "\t",
    BSgenomeName, 
    gRNA.file,
    outputDir,
    keep.R1only = TRUE, keep.R2only = TRUE,
    paired.direction = "opposite.strand",
    max.paired.distance = 1000, min.mapping.quality = 30,
    max.R1.len = 130, max.R2.len = 130,
    apply.both.max.len = FALSE, same.chromosome = TRUE,
    distance.inter.chrom = -1, min.R1.mapped = 30,
    min.R2.mapped = 30, apply.both.min.mapped = FALSE,
    max.duplicate.distance = 0,
    umi.plus.R1start.unique = TRUE, umi.plus.R2start.unique = TRUE,
    window.size = 20L, step = 20L, bg.window.size = 5000L,
    min.reads = 5L, min.SNratio = 2, maxP = 0.05,
    stats = c("poisson", "nbinom"),
    p.adjust.methods =
    c( "none", "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr"),
    distance.threshold = 40,
    plus.strand.start.gt.minus.strand.end = TRUE,
    gRNA.format = "fasta",
    overlap.gRNA.positions = c(17,18),
    upstream = 50, downstream = 50, PAM.size = 3, gRNA.size = 20,
    PAM = "NGG", PAM.pattern = "(NAG|NGG|NGA)$", max.mismatch = 6,
    allowed.mismatch.PAM = 2, overwrite = TRUE,
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,0.804, 0.685, 0.583))
}
\arguments{
  \item{alignment.inputfile}{
The alignment file. Currently supports bed output file with CIGAR information.
    Suggest run the workflow binReads.sh, which sequentially runs barcode
    binning, adaptor removal, alignment to genome, alignment quality filtering,
    and bed file conversion. Please download the workflow function and its
    helper scripts at http://mccb.umassmed.edu/GUIDE-seq/binReads/
}
  \item{umi.inputfile}{
A text file containing at least two columns, one is the read identifier and
    the other is the UMI or UMI plus the first few bases of R1 reads. Suggest
    use getUMI.sh to generate this file. Please download the script and its
    helper scripts at http://mccb.umassmed.edu/GUIDE-seq/getUMI/
}
  \item{alignment.format}{
The format of the alignment input file. Default bed file format. Currently only
    bed file format is supported, which is generated from binReads.sh
}
  \item{umi.header}{
Indicates whether the umi input file contains a header line or not. Default to
    FALSE
}
  \item{read.ID.col}{
The index of the column containing the read identifier in the umi input file,
    default to 1
}
  \item{umi.col}{
The index of the column containing the umi or umi plus the first few bases of
    sequence from the R1 reads, default to 2
}
  \item{umi.sep}{
column separator in the umi input file, default to tab
}
  \item{BSgenomeName}{
BSgenome object. Please refer to available.genomes in BSgenome package. For
    example, BSgenome.Hsapiens.UCSC.hg19 for hg19,
    BSgenome.Mmusculus.UCSC.mm10 for mm10,
    BSgenome.Celegans.UCSC.ce6 for ce6,
    BSgenome.Rnorvegicus.UCSC.rn5 for rn5,
    BSgenome.Drerio.UCSC.danRer7 for Zv9, and
    BSgenome.Dmelanogaster.UCSC.dm3 for dm3
}
 \item{gRNA.file}{
gRNA input file path or a DNAStringSet object that contains gRNA plus PAM
    sequences used for genome editing
}
  \item{outputDir}{
the directory where the off target analysis and reports will be written to
}
  \item{keep.R1only}{
Specify whether to include alignment with only R1 without paired R2.
    Default TRUE
}
  \item{keep.R2only}{
Specify whether to include alignment with only R2 without paired R1.
    Default TRUE
}
  \item{paired.direction}{
Specify whether the R1 and R2 should be aligned to the same strand or
    opposite strand. Default opposite.strand
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
Specify whether to apply minimum mapped length requirement to both R1 and R2
    reads, default FALSE
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
  \item{window.size}{
window size to calculate coverage
}
  \item{step}{
step size to calculate coverage
}
  \item{bg.window.size}{
window size to calculate local background
}
  \item{min.reads}{
minimum number of reads to be considered as a peak
}
  \item{min.SNratio}{
minimum signal noise ratio, which is the coverage normalized by local
    background
}
  \item{maxP}{
Maximum p-value to be considered as significant
}
  \item{stats}{
Statistical test, default poisson
}
  \item{p.adjust.methods}{
Adjustment method for multiple comparisons, default none
}
 \item{distance.threshold}{
Specify the maximum gap allowed between the plus stranded and
    the negative stranded peak, default 40. Suggest set it to twice of
    window.size used for peak calling.
}
  \item{plus.strand.start.gt.minus.strand.end}{
Specify whether plus strand peak start greater than
    the paired negative strand peak end. Default to TRUE
}
  \item{gRNA.format}{
Format of the gRNA input file. Currently, fasta is supported
}
  \item{PAM.size}{
PAM length, default 3
}
  \item{gRNA.size}{
The size of the gRNA, default 20
}
  \item{PAM}{
PAM sequence after the gRNA, default NGG
}
\item{overlap.gRNA.positions}{
The required overlap positions of gRNA and restriction enzyme cut site,
    default 17 and 18 for SpCas9.
}
  \item{max.mismatch}{
Maximum mismatch allowed in off target search, default 6
}
  \item{PAM.pattern}{
Regular expression of protospacer-adjacent motif (PAM), default
    (NAG|NGG|NGA)$ for off target search
}
\item{allowed.mismatch.PAM}{
Number of degenerative bases in the PAM sequence, default to 2 for N[A|G]G PAM
}
  \item{upstream}{
upstream offset from the peak start to search for off targets, default 50
}
  \item{downstream}{
downstream offset from the peak end to search for off targets, default 50
}
  \item{weights}{
a numeric vector size of gRNA length, default c(0, 0, 0.014, 0, 0, 0.395,
    0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,
    0.804, 0.685, 0.583) for SPcas9 system, which is used in Hsu et al., 2013
    cited in the reference section. Please make sure that the number of
    elements in this vector is the same as the gRNA.size, e.g., pad 0s at the
    beginning of the vector.
}
\item{overwrite}{
  overwrite the existing files in the output directory or not, default FALSE
}
}

\value{
\item{cleavage.gr }{Cleavage sites with one site per UMI as GRanges with
metadata column total set to 1 for each range}
\item{unique.umi.plus.R2}{a data frame containing unique cleavage site from
R2 reads mapped to plus strand with
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
   UMI (unique molecular identifier (umi) or umi with the first few bases of
   R1 read)
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
   UMI (unique molecular identifier (umi) or umi with the first few bases of
   R1 read)
}
\item{all.umi}{a data frame containing all the mapped reads with the
following columns.
   readName (read ID),   chr.x (chromosome of readSide.x/R1 read),
   start.x (start of eadSide.x/R1 read),
   end.x (end of eadSide.x/R1 read), mapping.qual.x
   (mapping quality of readSide.x/R1 read),
   strand.x (strand of readSide.x/R1 read),
   cigar.x (CIGAR of readSide.x/R1 read) , readSide.x (1/R1),
   chr.y (chromosome of readSide.y/R2 read)
   start.y (start of readSide.y/R2 read),  end.y (end of readSide.y/R2 read),
   mapping.qual.y (mapping quality of readSide.y/R2 read),
   strand.y (strand of readSide.y/R2 read),
   cigar.y (CIGAR of readSide.y/R2 read), readSide.y (2/R2),
   R1.base.kept (retained R1 length), R2.base.kept (retained R2 length),
   distance (distance between mapped R1 and R2),
   UMI (unique molecular identifier (umi) or umi with the first few bases
   of R1 read)
}
}
\references{
Shengdar Q Tsai and J Keith Joung et al. GUIDE-seq enables genome-wide
    profiling of off-target cleavage by CRISPR-Cas nucleases. Nature
    Biotechnology 33, 187 to 197 (2015)
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
        library("BSgenome.Hsapiens.UCSC.hg19")
        umiFile <- system.file("extdata", "UMI-HEK293_site4_R1.txt",
           package = "GUIDEseq")
        alignFile <- system.file("extdata","bowtie2.HEK293_site4.sort.bed" ,
            package = "GUIDEseq")
        gRNA.file <- system.file("extdata","gRNA.fa", package = "GUIDEseq")
        guideSeqRes <- GUIDEseqAnalysis(
            alignment.inputfile = alignFile,
            umi.inputfile = umiFile, gRNA.file = gRNA.file,
            BSgenomeName = Hsapiens, min.reads = 80)
        names(guideSeqRes)
   }
}
\keyword{ misc }

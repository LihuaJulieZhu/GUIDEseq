\name{GUIDEseqAnalysis}
\alias{GUIDEseqAnalysis}
\title{
Analysis pipeline for GUIDE-seq dataset
}
\description{
A wrapper function that uses the UMI sequence plus the first few bases of each
    sequence from R1 reads to estimate the starting sequence library, piles up
    reads with a user defined window and step size, identify the insertion 
    sites (proxy of cleavage sites), merge insertion sites from plus strand 
    and minus strand, followed by off target analysis of extended regions 
    around the identified insertion sites.
}
\usage{
GUIDEseqAnalysis(alignment.inputfile, umi.inputfile,
    alignment.format = c("auto", "bam", "bed"), 
    umi.header = FALSE, read.ID.col = 1L,
    umi.col = 2L, umi.sep = "\t",
    BSgenomeName, 
    gRNA.file,
    outputDir,
    n.cores.max = 3L,
    keep.R1only = TRUE, keep.R2only = TRUE,
    concordant.strand = TRUE,
    max.paired.distance = 1000L, min.mapping.quality = 30L,
    max.R1.len = 130L, max.R2.len = 130L,
    apply.both.max.len = FALSE, same.chromosome = TRUE,
    distance.inter.chrom = -1L, min.R1.mapped = 20L,
    min.R2.mapped = 20L, apply.both.min.mapped = FALSE,
    max.duplicate.distance = 0L,
    umi.plus.R1start.unique = TRUE, umi.plus.R2start.unique = TRUE,
    window.size = 20L, step = 20L, bg.window.size = 5000L,
    min.reads = 5L, min.reads.per.lib = 1L, 
    min.peak.score.1strandOnly = 5L,
    min.SNratio = 2, maxP = 0.05,
    stats = c("poisson", "nbinom"),
    p.adjust.methods =
    c( "none", "BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr"),
    distance.threshold = 40L,
    max.overlap.plusSig.minusSig = 30L,
    plus.strand.start.gt.minus.strand.end = TRUE,
    keepPeaksInBothStrandsOnly = TRUE, 
    gRNA.format = "fasta",
    overlap.gRNA.positions = c(17,18),
    upstream = 50L, downstream = 50L, PAM.size = 3L, gRNA.size = 20L,
    PAM = "NGG", PAM.pattern = "(NAG|NGG|NGA)$", max.mismatch = 6L,
    allowed.mismatch.PAM = 2L, overwrite = TRUE,
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,0.804, 0.685, 0.583),
    orderOfftargetsBy = c("peak_score", "predicted_cleavage_score", "n.mismatch"),
    descending = TRUE,
    keepTopOfftargetsOnly = TRUE,
    keepTopOfftargetsBy = c("predicted_cleavage_score", "n.mismatch"),
    scoring.method = c("Hsu-Zhang", "CFDscore"),
        subPAM.activity = hash( AA =0,
          AC =   0,
          AG = 0.259259259,
          AT = 0,
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
     txdb,
     orgAnn
)
}
\arguments{
  \item{alignment.inputfile}{
The alignment file. Currently supports bam and bed output file with CIGAR information.
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
\item{n.cores.max}{
Indicating maximum number of cores to use in multi core mode, 
    i.e., parallel processing, default 3. Please set it to 1 to disable
    multicore processing for small dataset.
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
\item{min.reads.per.lib}{
minimum number of reads in each library (usually two libraries)
to be considered as a peak
}
\item{min.peak.score.1strandOnly}{
Specify the minimum number of reads for a one-strand only peak to be included in the output.
Applicable when set keepPeaksInBothStrandsOnly to FALSE
}
  \item{min.SNratio}{
Specify the minimum signal noise ratio to be called as peaks, 
which is the coverage normalized by local background.
}
  \item{maxP}{
Specify the maximum p-value to be considered as significant
}
  \item{stats}{
Statistical test, currently only poisson is implemented
}
  \item{p.adjust.methods}{
Adjustment method for multiple comparisons, default none
}
 \item{distance.threshold}{
Specify the maximum gap allowed between the plus strand and
    the negative strand peak, default 40. Suggest set it to twice of
    window.size used for peak calling.
}
\item{max.overlap.plusSig.minusSig}{
Specify the cushion distance to allow sequence error and inprecise integration
Default to 30 to allow at most 10 (30-window.size 20) bp (half window) of 
minus-strand peaks on the right side of plus-strand peaks. Only applicable
if plus.strand.start.gt.minus.strand.end is set to TRUE.
}
  \item{plus.strand.start.gt.minus.strand.end}{
Specify whether plus strand peak start greater than
    the paired negative strand peak end. Default to TRUE
}
  \item{keepPeaksInBothStrandsOnly}{
Indicate whether only keep peaks present in both strands as specified
by plus.strand.start.gt.minus.strand.end, max.overlap.plusSig.minusSig and 
distance.threshold.
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
\item{overwrite}{
  overwrite the existing files in the output directory or not, default FALSE
}
  \item{weights}{
a numeric vector size of gRNA length, default c(0, 0, 0.014, 0, 0, 0.395,
    0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,
    0.804, 0.685, 0.583) for SPcas9 system, which is used in Hsu et al., 2013
    cited in the reference section. Please make sure that the number of
    elements in this vector is the same as the gRNA.size, e.g., pad 0s at the
    beginning of the vector.
}
\item{orderOfftargetsBy}{
 Criteria to order the offtargets, which 
 works together with the descending parameter 
}
\item{descending}{
 Indicate the output order of the offtargets, i.e., in the descending or ascending order. 
}
 \item{keepTopOfftargetsOnly}{
  Output all offtargets or the top offtarget using the keepOfftargetsBy criteria,
  default to the top offtarget
}
 \item{keepTopOfftargetsBy}{
  Output the top offtarget for each called peak using the keepTopOfftargetsBy criteria,
  If set to predicted_cleavage_score, then the offtargets with the highest predicted
  cleavage score will be retained 
  If set to n.mismatch, then the offtarget with the lowest number of mismatch to
  the target sequence will be retained 
}
 \item{scoring.method}{
   Indicates which method to use for offtarget cleavage rate estimation, currently
two methods are supported,  Hsu-Zhang and CFDscore
}
\item{subPAM.activity}{
   Applicable only when scoring.method is set to CFDscore
   A hash to represent the cleavage rate for each alternative sub PAM sequence relative
to preferred PAM sequence
}
\item{subPAM.position}{
   Applicable only when scoring.method is set to CFDscore
The start and end positions of the sub PAM. Default to 22 and 23 for SP with 20bp gRNA and NGG as preferred PAM
}
\item{PAM.location}{
PAM location relative to gRNA. For example, default to 3prime for spCas9 PAM. 
Please set to 5prime for cpf1 PAM since it's PAM is located on the 5 prime end
}
\item{mismatch.activity.file}{
Applicable only when scoring.method is set to CFDscore
A comma separated (csv) file containing the cleavage rates for all possible types of
single nucleotide mismatche at each position of the gRNA. By default, using the supplemental Table 19 from Doench et al., Nature Biotechnology 2016
}
\item{txdb}{
TxDb object, for creating and using TxDb object, please refer 
to GenomicFeatures package. For a list of existing TxDb object, please 
search for annotation package starting with Txdb at 
http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData, 
such as TxDb.Rnorvegicus.UCSC.rn5.refGene for rat, 
TxDb.Mmusculus.UCSC.mm10.knownGene for mouse, 
TxDb.Hsapiens.UCSC.hg19.knownGene for human, 
TxDb.Dmelanogaster.UCSC.dm3.ensGene for Drosophila and 
TxDb.Celegans.UCSC.ce6.ensGene for C.elegans
}
\item{orgAnn}{
organism annotation mapping such as org.Hs.egSYMBOL in org.Hs.eg.db package for human
}
}

\value{   
  \item{offTargets}{
a data frame, containing all input peaks
    with potential gRNA binding sites, mismatch number and positions, alignment
    to the input gRNA and predicted cleavage score.}
\item{merged.peaks}{merged peaks as GRanges
 with count (peak height), bg (local background), 
       SNratio (signal noise ratio), p-value, and option adjusted p-value }
\item{peaks }{GRanges with count (peak height), bg (local background), 
       SNratio (signal noise ratio), p-value, and option adjusted p-value }
\item{uniqueCleavages}{Cleavage sites with one site per UMI as GRanges with 
metadata column total set to 1 for each range}
\item{read.summary}{One table per input mapping file that contains the number of reads for each chromosome location}
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
if(!interactive())
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
            orderOfftargetsBy = "peak_score",
            descending = TRUE,
            keepTopOfftargetsBy = "predicted_cleavage_score",
            scoring.method = "CFDscore",
            BSgenomeName = Hsapiens, min.reads = 80, n.cores.max = 1)
        names(guideSeqRes)
   }
}
\keyword{ misc }

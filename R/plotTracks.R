#' Plot offtargets along all chromosomes with one track per chromosome
#'
#' @param offTargetFile The file path containing off-targets generated
#'  from GUIDEseqAnalysis
#' @param sep The separator in the file, default to tab-delimited
#' @param header Indicates whether the input file contains a header,
#' default to TRUE
#' @param gRNA.size The size of the gRNA, default 20
#' @param PAM.size PAM length, default 3
#' @param chromosome.order The chromosome order to plot from top to bottom
#' @param xlab The x-asix label, default to Chromosome Size (bp)
#' @param ylab The y-asix label, default to Peak Score. Change it to
#' be consistent with the score.col
#' @param score.col The column used as y values in the plot. Available choices
#' are peak_score, total.match, gRNA.match, total.mismatch.bulge,
#' gRNA.mismatch.bulge, and predicted_cleavage_score.
#' @param transformation Indicates whether plot the y-value in log10 scale or
#' in the original scale. When scale.col is set to total.match, gRNA.match,
#' total.mismatch.bulge, and gRNA.mismatch.bulge, the data will be plotted
#' in the original scale.
#' @param title The figure title, default to none.
#' @param axis.title.size The font size for the axis labels, default to 12
#' @param axis.label.size The font size for the tick labels, default to 8
#' @param strip.text.y.size The font size for the strip labels, default to 9
#' @param off.target.line.size The line size to depict the off-targets, default
#' to 0.6
#' @param on.target.line.size The line size to depict the on-targets, default
#' to 1
#' @param on.target.score The score for the on-target, default to 1 for CFD
#' scoring system. This is the maximum score in the chosen scoring system.
#' Change it accordingly if different off-target scoring system is used.
#' @param on.target.color The line color to depict the on-targets, default
#' to red
#' @param off.target.color The line color to depict the off-targets, default
#' to black
#' @param strip.text.y.angle The angel for the y strip text, default to 0.
#' Set it to 45 if angled representation is desired
#' @param scale.grid Used to set the scales in facet_grid, default to free_x,
#'  meaning that scales vary across different x-axis, but fixed in y-axis.
#'  Other options are fixed, free, and free_y meaning that scales shared
#'  across all facets, vary across both x- and y- axises, and vary across 
#'  y-axis only, respectively.
#'  For details, please type ?ggplot2::facet_grid
#'
#' @return a ggplot object

#' @importFrom ggplot2 ggplot aes theme scale_y_continuous
#' @importFrom ggplot2 element_blank element_text geom_rect
#' @importFrom ggplot2 theme_classic facet_grid xlab ylab ggtitle
#' @importFrom rlang sym
#' @importFrom dplyr '%>%'
#'
#' @export plotTracks
#' @author Lihua Julie Zhu
#'
#' @examples
#' if (interactive())
#' {
#'    offTargetFilePath <- system.file("extdata/forVisualization",
#'       "offTargetsInPeakRegions.xls",
#'        package = "GUIDEseq")
#'   fig1 <- plotTracks(offTargetFile = offTargetFilePath)
#'   fig1
#'   fig2 <- plotTracks(offTargetFile = offTargetFilePath,
#'     score.col = "total.mismatch.bulge",
#'     ylab = "Total Number of Mismatches and Bulges")
#'   fig2
#'   fig3 <- plotTracks(offTargetFile = offTargetFilePath,
#'      score.col = "total.match",
#'      ylab = "Total Number of Matches")
#'   fig3
#'   fig4 <- plotTracks(offTargetFile = offTargetFilePath,
#'       score.col = "gRNA.match",
#'       ylab = "Number of Matches in gRNA")
#'   fig4
#'   fig5 <- plotTracks(offTargetFile = offTargetFilePath,
#'       score.col = "gRNA.mismatch.bulge",
#'       ylab = "Number of Mismatches and Bulges in gRNA")
#'   fig5
#'   fig6 <- plotTracks(offTargetFile = offTargetFilePath,
#'      score.col = "predicted_cleavage_score",
#'      ylab = "CFD Score",
#'      scale.grid = "fixed",
#'      transformation = "none")
#'  fig6
#'  
#' }

plotTracks <- function(offTargetFile, sep ="\t",
                       header = TRUE,
                       gRNA.size = 20L,
                       PAM.size = 3L,
                       chromosome.order =
                         paste0("chr", c(1:22, "X", "Y")),
                       xlab  = "Chromosome Size (bp)",
                       ylab  = "Peak Score",
                       score.col = c("peak_score",
                                     "total.match",
                                     "gRNA.match",
                                     "total.mismatch.bulge",
                                     "gRNA.mismatch.bulge",
                                     "predicted_cleavage_score"
                                     ),
                       transformation = c("log10", "none"),
                       title = "",
                       axis.title.size = 12,
                       axis.label.size = 8,
                       strip.text.y.size = 9,
                       off.target.line.size = 0.6,
                       on.target.line.size = 1,
                       on.target.score = 1,
                       on.target.color = "red",
                       off.target.color = "black",
                       strip.text.y.angle = 0,
                       scale.grid = c( "free_x", "fixed", "free", "free_y")
                  )
{
   if(missing(offTargetFile) || !file.exists(offTargetFile))
      stop("Please specify the offTargetFile as a valid file path!")
   scale.grid <- match.arg(scale.grid)
   transformation <- match.arg(transformation)
   score.col <- match.arg(score.col)
   x <- read.table(offTargetFile,
                  sep = sep,
                  header = header,
                  stringsAsFactors = FALSE)
   x$chromosome <- factor(x$chromosome,
                          levels = chromosome.order)
   if(score.col == "peak_score" && min(x$peak_score)  < 2) # for log transformation
      x$peak_score <- x$peak_score  + 1
   else if (score.col == "predicted_cleavage_score" &&
            min(x$predicted_cleavage_score) == 0)
      x$predicted_cleavage_score[x$predicted_cleavage_score == 0] <-
         min(x$predicted_cleavage_score[x$predicted_cleavage_score > 0]) * 0.1

   # x1 <- makeGRangesFromDataFrame(x,
   #                                keep.extra.columns=TRUE,
   #                                ignore.strand=FALSE,
   #                                seqinfo=NULL,
   #                                seqnames.field="chromosome",
   #                                start.field = "offTarget_Start",
   #                                end.field = "offTarget_End",
   #                                strand.field = "offTargetStrand",
   #                                starts.in.df.are.0based=FALSE)

   if (score.col == "gRNA.mismatch.bulge")
      x <- x %>% mutate(gRNA.mismatch.bulge = 
                           total.mismatch.bulge - n.PAM.mismatch
                           )
   else if (score.col == "gRNA.match")
      x <- x %>% mutate(gRNA.match = gRNA.size -
                           (total.mismatch.bulge - n.PAM.mismatch)
      )
   else if (score.col == "total.match")
      x <- x %>% mutate(total.match = gRNA.size + PAM.size -
                           total.mismatch.bulge
      )
   
   score.col <- sym(score.col)
   p1 <- ggplot(x,
          aes(offTarget_Start, !!score.col)) +
   geom_rect(aes(xmin = offTarget_Start, ymin = 0,
                    xmax = offTarget_End,
                    ymax= !!score.col), color = off.target.color ,
                size = off.target.line.size) +
   geom_rect(data = x[x$predicted_cleavage_score == on.target.score,],
                         aes(xmin = offTarget_Start, ymin = 0,
                    xmax = offTarget_End,
                    ymax = !!score.col), color = on.target.color ,
                size = on.target.line.size
                )

    if (!as.character(score.col) %in% c("total.match",
        "gRNA.match",
        "total.mismatch.bulge",
        "gRNA.mismatch.bulge")  &&
        transformation == "log10")
        p1 <- p1 + scale_y_continuous(trans='log10',
                                      n.breaks = 4)
    else
        p1 <- p1 + scale_y_continuous(n.breaks = 4)

   if (length(unique(x$chromosome)) > 1)
      p1 <- p1 + facet_grid(chromosome ~.,
                scales = scale.grid )
   p1 <- p1 +
     theme_classic() +
      xlab(xlab) +
      ylab(ylab) +
      ggtitle(title) +
     theme(axis.title =element_text(size=axis.title.size,
                             face = "bold"),
           axis.text=element_text(size=axis.label.size),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           strip.text.y=element_text(angle=strip.text.y.angle,
                                     size = strip.text.y.size ),
           legend.position= "none"
        )
   p1
}

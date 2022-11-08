#' Plot offtargets as manhantann plots or along all chromosomes with one track per chromosome
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
#' are peak_score, n.distinct.UMIs, total.match, gRNA.match, total.mismatch.bulge,
#' gRNA.mismatch.bulge, and predicted_cleavage_score.
#' @param transformation Indicates whether plot the y-value in log10 scale or
#' in the original scale. When scale.col is set to total.match, gRNA.match,
#' total.mismatch.bulge, and gRNA.mismatch.bulge, transformation will
#' not be applied and the data will be plotted
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
#' @param  plot.type Plot type as tracks by individual chromosome or
#'  manhattan plot with all chromosome in one plot
#' @param family font family, default to sans (Arial). Other options are 
#' serif (Times New Roman) and mono (Courier). It is possible to use custom 
#' fonts with the extrafont package with the following commands
#' install.packages("extrafont")
#' library(extrafont)
#' font_import()
#' loadfonts(device = "postscript")
#' @param x.exp For transforming the x-axis to allow sufficient spaces
#' between small chromsoms default to 1.5
#' @param plot.zero.logscale Specifying "none" to filter out score.col with
#' zeros when plotting in log10 scale. Specify a very small numeric number
#' if you intend to show the zeros in log scale in the figure. If users
#' specify a number that's bigger than any positive score, plot.zero.logscale
#' will be set to the minimum positive score divided by 10.
#' @return a ggplot object
#'
#' @importFrom ggplot2 ggplot aes theme scale_y_continuous scale_color_manual scale_x_continuous
#' @importFrom ggplot2 element_blank element_text geom_rect element_line
#' @importFrom ggplot2 theme_classic facet_grid xlab ylab ggtitle element_blank
#' @importFrom rlang sym
#' @importFrom dplyr '%>%' summarise mutate left_join arrange select

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
#'  ## manhattan plot
#'   fig <- plotTracks(offTargetFile = offTargetFilePath,
#'         score.col = "total.mismatch.bulge", axis.title.size =9,
#'         plot.type =  "manhattan",
#'         ylab = "Number of Mismatches and Bulges in gRNA Plus PAM")
#'    fig
#'   fig <- plotTracks(offTargetFile = offTargetFilePath,
#'        score.col = "total.match", axis.title.size =9,
#'        plot.type =  "manhattan",
#'        ylab = "Number of Matches in gRNA Plus PAM")
#'    fig
#' fig <- plotTracks(offTargetFile = offTargetFilePath,
#'                  score.col = "gRNA.match",axis.title.size =9,
#'                  plot.type =  "manhattan",
#'                  ylab = "Number of Matches in gRNA")
#' fig
#' fig <- plotTracks(offTargetFile = offTargetFilePath,
#'                  score.col = "gRNA.mismatch.bulge", axis.title.size =9, 
#'                  plot.type =  "manhattan",
#'                  ylab = "Number of Mismatches and Bulges in gRNA")
#'   fig
#'   
#'   plotTracks(offTargetFile = offTargetFilePath,
#'       #'score.col = "predicted_cleavage_score",
#'       axis.title.size =9, family = "serif", plot.zero.logscale = 1e-6,
#'       plot.type =  "manhattan", transformation = "log10",
#'       ylab = "CFD Score")
#'       
#'   plotTracks(offTargetFile = offTargetFilePath,
#'        score.col = "peak_score",
#'        axis.title.size =9, 
#'        plot.type =  "manhattan",
#'        ylab = "Number of Insertion Events")
#'        
#'   plotTracks(offTargetFile = offTargetFilePath,
#'        score.col = "n.distinct.UMIs",
#'        axis.title.size =9, 
#'        plot.type =  "manhattan",
#'        ylab = "Number of Insertion Events")
#' }

plotTracks <- function(offTargetFile, sep ="\t",
                       header = TRUE,
                       gRNA.size = 20L,
                       PAM.size = 3L,
                       cleavage.position = 19L,
                       chromosome.order =
                         paste0("chr", c(1:22, "X", "Y")),
                       xlab  = "Chromosome Size (bp)",
                       ylab  = "Peak Score",
                       score.col = c("peak_score",
                                     "n.distinct.UMIs",
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
                       scale.grid = c( "free_x", "fixed", "free", "free_y"),
                       plot.type = c("manhattan", "tracks"),
                       family = "serif", x.exp = 1.5,
                       plot.zero.logscale = 1e-8
                  )
{
   if(missing(offTargetFile) || !file.exists(offTargetFile))
      stop("Please specify the offTargetFile as a valid file path!")
   scale.grid <- match.arg(scale.grid)
   transformation <- match.arg(transformation)
   score.col <- match.arg(score.col)
   
   score.col <- sym(score.col)
   plot.type <- match.arg(plot.type)
   x <- read.table(offTargetFile,
                  sep = sep,
                  header = header,
                  stringsAsFactors = FALSE)
   x$chromosome <- factor(x$chromosome,
                          levels = chromosome.order)
   
   # x1 <- makeGRangesFromDataFrame(x,
   #                                keep.extra.columns=TRUE,
   #                                ignore.strand=FALSE,
   #                                seqinfo=NULL,
   #                                seqnames.field="chromosome",
   #                                start.field = "offTarget_Start",
   #                                end.field = "offTarget_End",
   #                                strand.field = "offTargetStrand",
   #                                starts.in.df.are.0based=FALSE)

     x <- x %>% mutate(gRNA.mismatch.bulge = 
                           total.mismatch.bulge - n.PAM.mismatch,
                    gRNA.match = gRNA.size -
                           (total.mismatch.bulge - n.PAM.mismatch),
                    total.match = gRNA.size + PAM.size -
                           total.mismatch.bulge,
                    cleavage.position = ifelse(offTargetStrand == "+", 
                                               offTarget_Start + 
                                                 cleavage.position - 1,
                                               offTarget_End - 
                                                 cleavage.position + 1)
                    )
     
    if (plot.type == "manhattan") {
       x <- x %>%  group_by(chromosome) %>% 
          summarise(chr.max = max(cleavage.position)) %>% 
          mutate(chr.offset = cumsum(as.numeric(chr.max))-chr.max) %>%
          select(chromosome, chr.offset) %>%
          left_join(x, ., by=c("chromosome"="chromosome")) %>%
          arrange(chromosome, cleavage.position) %>%
          mutate(cum.cleavage.position = cleavage.position + chr.offset)
       
       ymin <- as.numeric(x %>% summarize(min(!!score.col))) - 2
       ymax <- as.numeric(x %>% summarize(max(!!score.col))) + 1
       xmax <- max(x$cum.cleavage.position)
       xmin <- min(x$cum.cleavage.position) - 1
       
       x$cum.cleavage.position <- 
         ((x$cum.cleavage.position - xmin)/xmin)^x.exp

       xaxis.lab.pos = x %>% group_by(chromosome) %>% 
             summarize(chr.center=( max(cum.cleavage.position) +
                           min(cum.cleavage.position) ) / 2 )

       xmin <- min(x$cum.cleavage.position) 
       if ( plot.zero.logscale == "none" && transformation == "log10")
           x <- x %>% filter(!!score.col <= 0)
       else if (transformation == "log10" && is.numeric(plot.zero.logscale)) {
           plot.zero.logscale <- min(plot.zero.logscale, 
                                     min(x[x[, as.character(score.col)] > 0, 
                                           as.character(score.col)]) / 10)
           x[x[, as.character(score.col)] <= 0, as.character(score.col)] <-
               plot.zero.logscale
       }
       # signed_log10 <- scales::trans_new("signed_log10",
       #                            transform=function(x) 
       #                               ifelse(x == 0, 0, sign(x)*log10(abs(x))),
       #                            inverse=function(x) 
       #                               ifelse(x == 0, 0, sign(x) * abs(x)^10))
       # 
       p1 <- ggplot(x, aes(x= cum.cleavage.position,
                           y = !!score.col)) +
             geom_point(aes(color=as.factor(chromosome)), 
                        size = off.target.line.size, shape = 6) +
             geom_point(data = x[x$predicted_cleavage_score == 
                                   on.target.score,],
                   aes(x = cum.cleavage.position,
                       y = !!score.col), shape = 25, 
                       fill = on.target.color,
                       color = on.target.color,
                       size = on.target.line.size) +
            scale_color_manual(values = rep(c("grey", "lightblue",
                                              "purple", "orange"),
                                            ceiling(nrow(xaxis.lab.pos)/4))) +
            scale_x_continuous(label = xaxis.lab.pos$chromosome, 
                                breaks= xaxis.lab.pos$chr.center) +
                              #  guide = guide_axis(n.dodge=n.dodge)) +
            theme_classic() 
       
      if (as.character(score.col) %in% c("total.match",
                                          "gRNA.match",
                                          "total.mismatch.bulge",
                                          "gRNA.mismatch.bulge"))
          p1 <- p1 + 
             scale_y_continuous(breaks = ymin:ymax) 
        else if (as.character(score.col) == "predicted_cleavage_score" && 
                transformation == "log10" && is.numeric(plot.zero.logscale)) {
            if (-log10(plot.zero.logscale) %% 2 == 0)
                breaks = 0.1^seq(0, -log10(plot.zero.logscale),
                                       by = 2)
            else
                breaks = 0.1^seq(0, -log10(plot.zero.logscale),
                               by = 1)
            labels <- breaks
            labels[length(labels)] = 0   
            p1 <- p1 + scale_y_continuous(trans='log10',
                                       breaks = breaks, labels = 
                                         format(labels, scientific = TRUE)) +
                theme(axis.line.y = element_blank()) +
                annotate(geom = "segment", x = -Inf, xend = -Inf,
                        # x = layer_scales(p1)$x$range$range[1], 
                        # xend = layer_scales(p1)$x$range$range[1], 
                         y = labels[length(labels) - 1], 
                         yend = Inf) +
                        # yend = 10^layer_scales(p1)$y$range$range[2]) +
                annotate(geom = "segment", x = -Inf, xend = -Inf,
                        # x = layer_scales(p1)$x$range$range[1], 
                        # xend = layer_scales(p1)$x$range$range[1], 
                         y =  labels[length(labels)], 
                         yend = labels[length(labels) - 1],
                    linetype = "dashed", color = "black")
        }
        else  if (transformation == "log10")
              p1 <- p1 + scale_y_continuous(trans='log10',n.breaks =6) 
            
      p1 <- p1 +
        #theme_bw() + 
        xlab(label = "") + 
        ylab(label = ylab) +
        theme( 
          legend.position="none",
          panel.border = element_blank(),
          axis.line = element_line(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title =element_text(size=axis.title.size,
                                   face = "bold" ,
                                   family = family),
                  axis.text=element_text(size=axis.label.size, 
                                         face = "bold" ,
                                         family = family),
                  axis.text.x=element_text(angle= 90, hjust= 1,
                                           vjust = 0.5,
                                           color = 
                                               rep(c("grey", "lightblue",
                                                      "purple", "orange"),
                                                    ceiling(nrow(xaxis.lab.pos)/4)))
              )
      if (as.character(score.col) %in% c("peak_score", "n.distinct.UMIs"))
        p1 <- p1 + theme(axis.text.y=element_text(angle=15,
                                                  hjust= 0.5,
                                 vjust = 0.5))
   }
   else {  ## plot.type is not manhattan
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
    p1 <- p1 + theme_bw() +
      xlab(xlab) +
      ylab(ylab) +
      ggtitle(title) +
      theme(axis.title =element_text(size=axis.title.size,
                             face = "bold", family = family),
           axis.text=element_text(size=axis.label.size, face = "bold",
                                  family = family),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           strip.text.y=element_text(angle=strip.text.y.angle,
                                     size = strip.text.y.size ),
           legend.position= "none"
        )
    } ### if plot.type is not manhattan
   p1
}

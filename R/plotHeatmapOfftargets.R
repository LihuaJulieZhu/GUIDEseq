#' Plot offtargets from multiple samples as heatmap
#'
#' @param mergedOfftargets a data frame from running the combineOfftargets 
#' function
#' @param min.detection.rate minimum relative detection rate to be included
#' in the heatmap
#' @param font.size font size for x labels and numbers along the y-axis.
#' @param on.target.predicted.score Default to 1 for the CFDscore
#' scoring method. Set it to 100 for the Hsu-Zhang
#' scoring method.
#' @param IR.normalization Default to sequence.depth which uses
#' the sequencing depth for each sample
#' in the input file to calculate the relative insertion rate (RIR). 
#' Other options are "on.target.score" and "sum.score" which use the 
#' on-target score for each sample and the sum of all on-target and 
#' off-target scores to calculate the RIR respectively.
#' The score can be either peak.score or n.distinct.UMIs as specified by
#' the parameter insertion.score.column
#' @param top.bottom.height.ratio the ratio of the height of top panel
#' vs that of the bottom panel.
#' @param dot.distance.breaks a numeric vector for specifying
#' the minimum number of rows in each panel to use the
#' the corresponding distance in dot.distance.scaling.factor
#' between consecutive dots along the y-axis.
#' In the default setting, dot.distance.breaks and
#' dot.distance.scaling.factor are
#' set to c(5, 10, 20, 40, 60) and c(0.4, 0.6, 0.8, 1.2, 2)
#' respectively, which means that if the number of rows in each panel
#' is greater than or equal to 60, 40-59, 20-39, 10-19, 5-9, and less
#' than 5,then the distance between consecutive dots will
#' be plotted 2, 1.2, 0.8, 0,6, 0.4, and 0.2 (half of 0.4) units away
#' in y-axis respectively.
#' @param dot.distance.scaling.factor a numeric vector for specifing
#' the distance between two consecutive dots.
#' See dot.distance.breaks for more information.
#' @param bottom.start.offset Default to 2, means that place the
#' top number in the bottom panel 2 units below the top border. Increase
#' the value will move the number away from the top border.
#' @param color.low The color used to represent the lowest indel rate,
#' default to white
#' @param color.high The color used to represent the highest indel rate
#' the intermediate indel rates will be colored using the color between
#' color.low and color.high. Default to blue.
#' @param sample.names Optional sample Names used to label the x-axis. If
#' not provided, x-axis will be labeled using the sample names provided
#' in the GUIDEseqAnalysis step.
#' @param insertion.score.column "n.distinct.UMI" or "peak_score" to be included on
#' the right side of the alignment as Insertion Events. Relative Insertion 
#' Rate (RIR) % is calculated as peak_score/n.distinct.UMI 
#' divided by ontarget peak_score/n.distinct.UMI. For example, RIR for ontarget
#' should be 100
#' @return a ggplot object
#' @importFrom ggplot2 ggplot aes annotate theme geom_tile geom_point
#' @importFrom ggplot2 element_blank element_text element_rect
#' @importFrom ggplot2 scale_fill_gradient theme_bw
#' @importFrom purrr map
#' @importFrom dplyr group_by arrange filter mutate '%>%'
#' @importFrom tidyr unnest
#' @import patchwork
#' @export plotAlignedOfftargets
#' @author Lihua Julie Zhu
#'
#' @examples
#' if (interactive())
#' {
#'   mergedOfftargets <- 
#'         read.table(system.file("extdata/forVisualization",
#'       "mergedOfftargets.txt",
#'        package = "GUIDEseq"),
#'                    sep ="\t", header = TRUE)
#'                    
#'  figs <- plotHeatmapOfftargets(mergedOfftargets,
#'                    min.detection.rate = 2.5,
#'                    IR.normalization = "on.target.score",
#'                    top.bottom.height.ratio = 12,
#'                    bottom.start.offset = 6,
#'                    dot.distance.scaling.factor = c(0.2,0.2,0.4,0.4, 0.4),
#'                    sample.names = c("Group1", "Group2"))
#'                    figs[[1]]/figs[[2]] +
#'  plot_layout(heights = unit(c(2,1),
#'                              c('null', 'null')))
#'                              
#' figs = plotHeatmapOfftargets(mergedOfftargets,
#'                  min.detection.rate = 1.2,
#'                  IR.normalization = "sum.score",
#'                  top.bottom.height.ratio = 12,
#'                  bottom.start.offset = 6,
#'                  dot.distance.scaling.factor = c(0.2,0.2,0.4,0.4, 0.4),
#'                  sample.names = c("Group1", "Group2"))
#'                  figs[[1]]/figs[[2]] +
#'                  plot_layout(heights = unit(c(2,1),
#'                   c('null', 'null')))
#'  figs <- plotHeatmapOfftargets(mergedOfftargets,
#'     min.detection.rate = 0.2,
#'     IR.normalization = "sequence.depth",
#'     top.bottom.height.ratio = 12,
#'     bottom.start.offset = 6,
#'     dot.distance.scaling.factor = c(0.2,0.2,0.2,0.2, 0.2),
#'     sample.names = c("Group1", "Group2"))
#' figs[[1]]/figs[[2]] +
#'     plot_layout(heights = unit(c(2,1),
#'     c('null', 'null')))
#' figs = plotHeatmapOfftargets(mergedOfftargets,
#'     min.detection.rate = 3,
#'     IR.normalization = "none",
#'     top.bottom.height.ratio = 12,
#'     bottom.start.offset = 6,
#'     dot.distance.scaling.factor = c(0.2,0.2,0.7,0.7, 0.7),
#'     sample.names = c("Group1", "Group2"))
#'     figs[[1]]/figs[[2]] 
#' plot_layout(heights = unit(c(2,1),
#'                 c('null', 'null')))
#' }

plotHeatmapOfftargets <- function(mergedOfftargets,
                           min.detection.rate = 0.1,
                           font.size = 12,
                           on.target.predicted.score = 1,
                           IR.normalization = c("sequence.depth", 
                                                "on.target.score", 
                                                "sum.score", "none"),
                           top.bottom.height.ratio = 3,
                           dot.distance.breaks =
                             c(5, 10, 20, 40, 60),
                           dot.distance.scaling.factor =
                             c(0.4, 0.6, 0.8, 1.2, 2),
                           bottom.start.offset = 8,
                           color.low = "white",
                           color.high = "blue",
                           sample.names,
                           insertion.score.column = c("n.distinct.UMIs",
                                                      "peak_score"
                           )
                        )
{
  if(missing(mergedOfftargets) ||
     class(mergedOfftargets) !=  "data.frame")
   stop("mergedOfftargets is required as a data frame!")

  IR.normalization <- match.arg(IR.normalization)
  
  insertion.score.column <- match.arg(insertion.score.column)
  
  score.col <- grep(insertion.score.column,
                          colnames(mergedOfftargets))
  if (length(score.col) == 0)
  {
    insertion.score.column <- "peak_score"
     score.col <- grep(insertion.score.column,
                      colnames(mergedOfftargets))
  }
    
  if (missing(sample.names))
    sample.names <- gsub(paste0(".", insertion.score.column), "",
                       colnames(mergedOfftargets)[score.col],
                       fixed = TRUE)
  
  insertion.score.column <- sym(insertion.score.column)
  
  on.target.IR <- mergedOfftargets %>% 
    filter(total.mismatch.bulge == 0)
  on.target.IR <-
    rowMaxs(t(as.matrix((on.target.IR[, score.col]))), na.rm = TRUE)
  
  x <- mergedOfftargets[,score.col]

  for (i in 1:length(score.col))
  {
    x[is.na(x[,i]), i] <- 0
  }

  ontarget.dot.expand <- 0.55

  ontarget.row <- which(mergedOfftargets$predicted_cleavage_score ==
                          on.target.predicted.score)
  ontarget = rep(0, nrow(x))
  ontarget[ontarget.row] <- 1
  ontarget <- unlist(map(ontarget, rep,
             length(sample.names)))


  if(IR.normalization == "sequence.depth")
  {
    seq.depth.col <- grep("sequence.depth",
                      colnames(mergedOfftargets))
    seq.depth <- rowMaxs(t(as.matrix(
       mergedOfftargets[,seq.depth.col])), na.rm = TRUE)
    x2 <- round(t(t(x)/seq.depth) * 100, digits = 3)
  }
  else if (IR.normalization == "sum.score")
  {
      x2 <- round(t(t(x)/colSums(x)) * 100, digits = 3)
  }
  else if (IR.normalization == "on.target.score")
  {
      x2 <- round(t(t(x)/on.target.IR) * 100, digits = 3)
  }
  else
  {
      x2 <- round(log10(x+1), digits = 3)
  }

  y.labs = unlist(map(1:nrow(x), rep,
                          length(sample.names)))
  x0 <- data.frame(Offtargets = y.labs,
                   Samples = rep(sample.names, nrow(x)),
                   IR = as.vector(t(x2)), Ontarget = ontarget) %>%
    group_by(Offtargets) %>%
    mutate(IR.all = sum(IR), IR.max = max(IR)) %>%
    arrange(desc(Ontarget), desc(IR.max))

  x0$Offtargets <- sort(y.labs, decreasing = TRUE)

  x0.top <- x0 %>%
    filter(max(IR) >= min.detection.rate)
  x0.bottom <- x0 %>%
    filter(max(IR) < min.detection.rate) %>%
    mutate(IR = 0)

  dot.distance.m <- data.frame(dot.distance.breaks =
                                 dot.distance.breaks,
                               dot.distance.scaling.factor =
                                 dot.distance.scaling.factor) %>%
    arrange(dot.distance.breaks)

  if (nrow(x0.bottom) > 0)
  {
    pb <- ggplot(x0.bottom,
           mapping = aes(x = Samples,
                         y = Offtargets)) +
         geom_tile(fill = "white") +
        # geom_rect(data = x0.bottom, fill = "white",
        #       mapping = aes(xmin = 1,
        #                     xmax = length(sample.names),
        #                     ymin = 1,
        #                     ymax = max(Offtargets))) +
    annotate("text", x = 0.42,
             y = max(x0.bottom$Offtargets) - bottom.start.offset,
             label = max(x0$Offtargets) - nrow(x0.bottom)/length(sample.names) + 1,
             hjust = "left",
             size = font.size/3,
             vjust = "bottom") +
    theme_bw() +
    theme(text=element_text(size=font.size, family="sans",
                            face = "bold"),
          axis.text.x=element_text(size = font.size,
                                   family="sans",
                                   face = "bold",
                                   hjust = 0.6),
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", size = 1),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank())

    if (nrow(x0.bottom)/length(sample.names) <= 5)
      dist.dot <- min(dot.distance.scaling.factor)/2 *
        top.bottom.height.ratio
    else if (nrow(x0.top) > 5)
      dist.dot <-  dot.distance.m[max(which(
        dot.distance.m[,1] < nrow(x0.top))),2] *
        nrow(x0.bottom) / nrow(x0.top) * top.bottom.height.ratio
    else
      dist.dot <- dot.distance.m[max(which(
          dot.distance.m[,1] < nrow(x0.bottom))),2]

    if (nrow(x0.bottom)/length(sample.names) > 1)
      pb <- pb + annotate("text", x = 0.42, y = 1,
             label = max(x0$Offtargets), hjust = "left",
             size = font.size/3, vjust = "bottom")
    if (nrow(x0.bottom)/length(sample.names) >= 5)
      pb <- pb +
          annotate("text", x = 0.44,
             y = mean(x0.bottom$Offtargets) - dist.dot,
             label = ".", hjust = "left", size = font.size/2,
             vjust = "bottom") +
      annotate("text", x = 0.44, y = mean(x0.bottom$Offtargets),
             label = ".", hjust = "left", size = font.size/2,
             vjust = "bottom")  +
      annotate("text", x = 0.44,
             y = mean(x0.bottom$Offtargets) + dist.dot,
             label = ".", hjust = "left", vjust = "bottom",
             size = font.size/2) 
  }
  if (nrow(x0.top) > 0)
  {
     pu <- ggplot(x0.top,
               mapping = aes(x = Samples,
                             y = Offtargets, fill = IR)) +
            geom_tile() +
            guides(fill=guide_legend(title=NULL)) +
            geom_point(aes(x = length(sample.names) +
                             ontarget.dot.expand,
                    y =  x0.top$Offtargets[x0.top$Ontarget == 1][1]),
                    size = 2,
                    colour="red",
                    show.legend = FALSE) +
            annotate("text", x = 0.44, y = max(x0.top$Offtargets),
               label = "1", hjust = "left", size = font.size/3,
               vjust = "bottom") +
     # annotate("text", x = length(sample.names) + 0.5,
     #         y = x0.top$Offtargets[x0.top$Ontarget == 1][1],
     #         label = ".", hjust = "left", vjust = "top",
     #         size = 16, color ="red")+
          theme_bw()  +
          scale_fill_gradient(low=color.low, high= color.high) 

     if (nrow(x0.top)/length(sample.names) <= 5)
       dist.dot <- min(dot.distance.scaling.factor)/2
     else
       dist.dot <- dot.distance.m[max(which(
         dot.distance.m[,1] < nrow(x0.top))),2]
    if (nrow(x0.top)/length(sample.names) > 1)
      pu <- pu + annotate("text", x = 0.42, y = min(x0.top$Offtargets),
              label = nrow(x0.top)/length(sample.names), hjust = "left",
              size = font.size/3, vjust = "bottom")
    if (nrow(x0.top)/length(sample.names) >= 5)
      pu <- pu +
          annotate("text", x = 0.44,
                   y = mean(x0.top$Offtargets) - dist.dot,
               label = ".", hjust = "left", size = font.size/2,
               vjust = "bottom") +
          annotate("text", x = 0.44, y = mean(x0.top$Offtargets),
               label = ".", hjust = "left", size = font.size/2,
               vjust = "bottom")  +
          annotate("text", x = 0.44,
                   y = mean(x0.top$Offtargets) + dist.dot,
               label = ".", hjust = "left", vjust = "bottom",
               size = font.size/2)
  }
  if (nrow(x0.top) > 0 && nrow(x0.bottom) >0 )
  {
    pu <- pu +
      theme(text=element_text(size = font.size, family="sans",
                              face = "bold"),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", size = 1),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_blank()
      )
    print(pu/pb +
      plot_layout(widths = c(2, 1),
                  heights = unit(c(top.bottom.height.ratio,1),
                                 c('null', 'null'))))
    list(pu, pb)
  }
  else if (nrow(x0.top) > 0)
  {
    pu <- pu +
    theme(text=element_text(size = font.size, family="sans",
                            face = "bold"),
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", size = 1),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank()
    )
    pu
  }
  else
    pb
}

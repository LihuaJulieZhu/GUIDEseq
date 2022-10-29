#' Plot offtargets aligned to the target sequence
#'
#' @param offTargetFile The path of the file offTargetsInPeakRegions.xls
#' that stores the offtargets to be plotted. This file is the
#' output file from the function GUIDEseqAnalysis.
#' @param sep Field delimiter for the file specified as
#' offTargetFile, default to tab dilimiter
#' @param header Indicates whether there is header in the file
#' specified as offTargetFile, default to TRUE
#' @param gRNA.size Size of the gRNA, default to 20 for SpCas9 system
#' @param input.DNA.bulge.symbol The symbol used to represent DNA bulges
#' in the file specified as offTargetFile, default to "^"
#' @param input.RNA.bulge.symbol The symbol used to represent RNA bulges
#' in the file specified as offTargetFile, default to "-"
#' @param input.match.symbol The symbol used to represent matched bases
#' in the file specified as offTargetFile, default to "."
#' @param plot.DNA.bulge.symbol The symbol used to represent DNA bulges
#' in the figure to be generated, default to "I"
#' @param plot.RNA.bulge.symbol The symbol used to represent RNA bulges
#' in the figure to be generated, default to " "
#' @param plot.match.symbol The symbol used to represent matched bases
#' in the figure to be generated, default to "-"
#' @param color.DNA.bulge The color used to represent DNA bulges
#' in the figure to be generated, default to "red"
#' @param size.symbol The size used to plot the bases, and the symbols
#' of DNA/RNA bulges, default to 3
#' @param color.values The color used to represent different bases, DNA
#' bulges, and RNA bulges.
#' @param PAM PAM sequence in the target site, please update it to
#' the exact PAM sequence in the input target site.
#' @param body.tile.height  Specifies the height of each plotting
#' tile around each base/symbol for offtargets, default to 2.5
#' @param header.tile.height Specifies the height of each plotting
#' tile around each base/symbol for the target sequence on the very
#' top, default to 4
#' @param hline.offset Specifies the offset from the top border to draw
#' the horizontal line below the gRNA sequence, default to 4. Increase it
#' to move the line down and decrease it to move the line up.
#' @param plot.top.n Optional. If not specified, all the offtargets
#' in the input file specified as offTargetFile will be included
#' in the plot. For samples with a very large number of offtargets,
#' users can select the top n offtargets to be included in the plot.
#' For example, set plot.top.n = 20 to include only top 20 offtargets
#' in the plot. Please note offtargets are ordered by the peak_score
#' from top to bottom. Insertion Rate (IR) % is calculated as peak_score divided
#' by sum(peak_score) of all offtargets in the offtarget file.
#'
#' @return a ggplot object
#' @importFrom ggplot2 ggplot aes geom_text annotate theme geom_tile
#' @importFrom ggplot2 element_blank element_text element_rect
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme_bw
#' @importFrom ggplot2 geom_hline geom_vline scale_fill_manual
#' @importFrom purrr map
#' @importFrom dplyr select arrange mutate top_n '%>%'
#' @importFrom tidyr unnest
#' @export plotAlignedOfftargets
#' @author Lihua Julie Zhu
#'
#' @examples
#' offTargetFilePath <- system.file("extdata/forVisualization",
#'  "offTargetsInPeakRegions.xls",
#'  package = "GUIDEseq")
#' fig1 <- plotAlignedOfftargets(offTargetFile = offTargetFilePath,
#'     plot.top.n = 20)
#' fig1
#'
plotAlignedOfftargets <- function(offTargetFile, sep ="\t",
                           header = TRUE, gRNA.size = 20L,
                           input.DNA.bulge.symbol = "^",
                           input.RNA.bulge.symbol = "-",
                           input.match.symbol = ".",
                           plot.DNA.bulge.symbol = "I",
                           plot.RNA.bulge.symbol = " ",
                           plot.match.symbol = "-",
                           color.DNA.bulge = "red",
                           size.symbol = 3,
                           color.values =
                             c("A" = "#B5D33D", #green
                               "T" = "#AE9CD6", #purple
                               "C" = "#6CA2EA", #blue
                               "G" = "#FED23F", #orange
                               " " = "gray",
                               "-" = "white"),
                           PAM = "GGG",
                           body.tile.height = 2.5,
                           header.tile.height = 4,
                           hline.offset = 4,
                           plot.top.n
                        )
{
  if(missing(offTargetFile) || !file.exists(offTargetFile))
   stop("Please specify the offTargetFile as a valid file path!")
  x <- read.table(offTargetFile,
                  sep = sep,
                  header = header,
                  stringsAsFactors = FALSE)
  x <- x %>% mutate(IR =
                      round(peak_score/sum(peak_score) * 100,
                            digits = 2)) %>%
                      arrange(IR)
  if (!missing(plot.top.n) && class(plot.top.n) == "numeric")
    x <- x %>% top_n(plot.top.n, IR)

  x <- rbind(x[x$total.mismatch.bulge > 0,], x[x$total.mismatch.bulge == 0,])
  ontarget <- paste0(substr(x$gRNAPlusPAM[1],
                                            1, gRNA.size),
                                     PAM)
  aligns <- x %>% select(guideAlignment2OffTarget) %>%
    as.vector %>% unlist

  PAM.char <- matrix(unlist(map(x$PAM.sequence, strsplit, "")),
                     nrow = nrow(x), byrow = TRUE)
  PAM.ontarget <- unlist(strsplit(PAM, ""))
  PAM.char[PAM.char[,1] == PAM.ontarget[1], 1] <- input.match.symbol
  PAM.char[PAM.char[,2] == PAM.ontarget[2], 2] <- input.match.symbol
  PAM.char[PAM.char[,3] == PAM.ontarget[3], 3] <- input.match.symbol

  aligns <- c(paste0(aligns,
                   paste0(PAM.char[,1],
                          PAM.char[,2],
                          PAM.char[,3])),
              ontarget)

  aligns.char <-  unlist(map(aligns, strsplit, ""))
  aligns.char <- aligns.char[aligns.char != input.DNA.bulge.symbol]

  aligns.char[aligns.char == input.RNA.bulge.symbol] <- plot.RNA.bulge.symbol
  aligns.char[aligns.char == input.match.symbol] <- plot.match.symbol


  x0 <- data.frame(x = rep(1:(gRNA.size + nchar(PAM)),
                           length(aligns) ),
              y = 2 * sort(rep(1:length(aligns),
                               gRNA.size + nchar(PAM))),
              aligns.char = aligns.char,
              h = c(rep(body.tile.height, nrow(x) * (gRNA.size + nchar(PAM))),
                rep(header.tile.height, gRNA.size + nchar(PAM))))
  x0.for.symbol <- x0 %>%
    mutate(y = y - 0.6)
  x1 <- data.frame(x = rep(gRNA.size + nchar(PAM) + 1.6,
                           nrow(x)),
                   y = 2 * 1:nrow(x) - 0.6,
                   IR = x$IR)
  x.DNA.bulge <- data.frame(x = x$pos.DNA.bulge,
                   y = 2 * 1:nrow(x))
  x.DNA.bulge <- x.DNA.bulge[!is.na(x.DNA.bulge[,1]) &
                               x.DNA.bulge[,1] != "",] %>%
    mutate(x = strsplit(as.character(x), ",")) %>%
    unnest(x) %>%
    mutate(x = as.numeric(x) - 0.5)

  max.y <- max(x0$y) + 2
  max.x <- max(x0$x)
  ggplot(x0, aes(x = x, y = y)) +
     geom_tile(aes(x, y, fill = aligns.char,
                   height=h)) +
     scale_fill_manual(values = color.values) +
     geom_text(data = x0.for.symbol, aes(x=x,
                               y=y, label = aligns.char),
               size = size.symbol) +
     geom_text(data = x1, aes(x, y, label = IR), hjust = "middle") +
     annotate("text", x = max(x1$x), y = max.y - 2,
              label = "IR (%)", hjust = "middle") +
     #annotate("text", x = x1$x, y = x1$y, label = x1$IR) +
     geom_text(data = x.DNA.bulge, aes(x=x,
                             y=y, 
                             label = plot.DNA.bulge.symbol),
                             size = size.symbol + 4,
              color = color.DNA.bulge) +
     geom_vline(xintercept = gRNA.size + 0.5) +
     geom_vline(xintercept = gRNA.size + nchar(PAM) + 0.5) +
     geom_hline(yintercept = max.y - hline.offset) +
     theme_bw() +
     scale_y_continuous(limits = c(0, max.y),
                        expand = c(0, 0)) +
     scale_x_continuous(limits = c(0.5, max(x1$x) + 1.4),
                       expand = c(0, 0)) +
     #coord_cartesian(xlim = c(0, max(x1$x) + 1), clip = "off") +
     theme(text=element_text(size=size.symbol, family="Arial",
                            face = "bold"),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", size = 2),
          legend.position= "none",
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank()
          )
}

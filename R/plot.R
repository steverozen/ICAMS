#' Plot \strong{one} spectrum or signature.
#'
#' Plot the spectrum of \strong{one} sample or plot \strong{one} signature. The
#' type of graph is based on one attribute("catalog.type") of the input catalog.
#' You can first use \code{\link{TransformCatalog}} to get different types of
#' catalog and then do the plotting.
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}} with attributes added.
#' See \code{\link{as.catalog}} for more details.
#'
#' @param plot.SBS12 Only meaningful for class \code{SBS192Catalog}; if \code{TRUE},
#' generate an abbreviated plot of only SBS without context, i.e. 
#' C>A, C>G, C>T, T>A, T>C, T>G each on transcribed and untranscribed strands,
#' rather than SBS in trinucleotide context, e.g.
#' ACA > AAA, ACA > AGA, ..., TCT > TAT, ...
#'
#' @param cex Has the usual meaning. Taken from \code{par("cex")} by default.
#'   Only implemented for SBS96Catalog, SBS192Catalog and DBS144Catalog.
#'
#' @param grid A logical value indicating whether to draw grid lines. Only
#'   implemented for SBS96Catalog.
#'
#' @param upper A logical value indicating whether to draw horizontal lines and
#'   the names of major mutation class on top of graph. Only implemented for
#'   SBS96Catalog.
#'
#' @param ylim Has the usual meaning. Only implemented for SBS96Catalog.
#'
#' @param xlabels A logical value indicating whether to draw x axis labels. Only
#'   implemented for \code{SBS96Catalog}. If \code{FALSE} then plot x axis tick marks;
#'   set \code{par(tck = 0)} to suppress.
#'
#' @import graphics
#'
#' @return A list whose first element is a logic value indicating whether the
#'   plot is successful. The second element is a numeric vector giving the
#'   coordinates of all the bar midpoints drawn, useful for adding to the
#'   graph(\strong{currently only implemented for SBS96Catalog}).
#'
#' @note The sizes of repeats involved in deletions range from 0 to 5+ in the
#'   mutational-spectra and signature catalog rownames, but for plotting and
#'   end-user documentation deletion repeat sizes range from 1 to 6+.
#'
#' @export
#'
#' @name PlotCatalog
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "strelka.regress.cat.sbs.96.csv",
#'                     package = "ICAMS")
#' catSBS96 <- ReadCatalog(file)
#' colnames(catSBS96) <- "sample"
#' PlotCatalog(catSBS96)
PlotCatalog <- function(catalog, plot.SBS12 = NULL, cex = NULL, 
                        grid = NULL , upper = NULL, xlabels = NULL,
                        ylim = NULL) {
  UseMethod(generic = "PlotCatalog")
}

#' Plot catalog to a PDF file.
#'
#' Plot catalog to a PDF file. The type of graph is based on one
#' attribute("catalog.type") of the input catalog. You can first use
#' \code{\link{TransformCatalog}} to get different types of catalog and then do
#' the plotting.
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}} with attributes added.
#' See \code{\link{as.catalog}} for more details.
#'
#' @param file The name of the PDF file to be produced.
#'
#' @param plot.SBS12 Only meaningful for class \code{SBS192Catalog}; if \code{TRUE},
#' generate an abbreviated plot of only SBS without context, i.e. 
#' C>A, C>G, C>T, T>A, T>C, T>G each on transcribed and untranscribed strands,
#' rather than SBS in trinucleotide context, e.g.
#' ACA > AAA, ACA > AGA, ..., TCT > TAT, ... There are 12 bars in the graph.
#'
#' @param cex Has the usual meaning. A default value has been used by the
#'   program internally. Only implemented for SBS96Catalog, SBS192Catalog and
#'   DBS144Catalog.
#'   
#' @param grid A logical value indicating whether to draw grid lines. Only
#'   implemented for SBS96Catalog.
#'
#' @param upper A logical value indicating whether to draw horizontal lines and
#'   the names of major mutation class on top of graph. Only implemented for
#'   SBS96Catalog.
#'
#' @param xlabels A logical value indicating whether to draw x axis labels. Only
#'   implemented for \code{SBS96Catalog}. If \code{FALSE} then plot x axis tick marks;
#'   set \code{par(tck = 0)} to suppress.
#'   
#' @param ylim Has the usual meaning. Only implemented for SBS96Catalog.
#'
#' @return \code{invisible(TRUE)}
#'   
#' @note The sizes of repeats involved in deletions range from 0 to 5+ in the
#'   mutational-spectra and signature catalog rownames, but for plotting and
#'   end-user documentation deletion repeat sizes range from 1 to 6+.
#'
#' @export
#'
#' @name PlotCatalogToPdf
#' 
#' @examples 
#' file <- system.file("extdata",
#'                     "strelka.regress.cat.sbs.96.csv",
#'                     package = "ICAMS")
#' catSBS96 <- ReadCatalog(file)
#' colnames(catSBS96) <- "sample"
#' PlotCatalogToPdf(catSBS96, file = file.path(tempdir(), "test.pdf"))
PlotCatalogToPdf <- 
  function(catalog, file, plot.SBS12 = NULL, 
           cex = NULL, grid = NULL, 
           upper = NULL, 
           xlabels = NULL, 
           ylim = NULL) {
  UseMethod(generic = "PlotCatalogToPdf")
}

###############################################################################
# Plotting functions for SBS96, SBS192 and SBS1536 catalog start here
###############################################################################

#' @export
PlotCatalog.SBS96Catalog <-
  function(catalog, plot.SBS12, cex = par("cex"), grid = TRUE,
           upper = TRUE, xlabels = TRUE, ylim = NULL) {
    # stopifnot(dim(catalog) == c(96, 1))
    stopifnot(rownames(catalog) == ICAMS::catalog.row.order$SBS96)

    class.col <- c("#0000ff",  
                   "#000000",  
                   "#ff4040",  
                   "#838383",  
                   "#40ff40",  
                   "#ff667f")  
    
    if (ncol(catalog) == 1) {
      cols <- rep(class.col, each = 16)
      to.plot <- catalog[ ,1]
    } else if (ncol(catalog) == 2) {
      cols <- c("red", "grey35")
      to.plot <- t(catalog) # [ , 1]
    } else {
      stop("Can only handle 1 or 2 column catalogs")
    }
    maj.class.names <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    num.classes <- 96

    catalog.type <- attributes(catalog)$catalog.type
    if (catalog.type == "density") {
      ylab <- "mut/million"
      to.plot <- 1e6 * to.plot
      ymax    <- 1e6 * max(rowSums(catalog))
    } else if (catalog.type == "counts") {
      # Set a minimum value for ymax to make the plot more informative
      ymax <- 4 * ceiling(max(max(rowSums(catalog)), 10) / 4)
      ylab <- "counts"
    } else if (catalog.type %in%
               c("counts.signature", "density.signature")) {
      ylab <- ifelse(catalog.type == "counts.signature",
                     "counts proportion", "density proportion")
      ymax <- max(rowSums(catalog))
    } else {
      stop("Programming error, illegal catalog type ", catalog.type)
    }
    if (is.null(ylim)) {
      ylim <- c(0, ymax)
    } else {
      ymax <- ylim[2]
    }
      
    bp <- barplot(to.plot, xaxt = "n", yaxt = "n", xaxs = "i",
                  xlim = c(-1, 230),
                  ylim = ylim, lwd = 3, space = 1.35, border = NA,
                  col = cols, ylab = ylab, cex.lab = cex * par("cex.lab"))

    # Draw the x axis
    segments(bp[1] - 0.5, 0, bp[num.classes] + 0.5, 0, 
             col = 'grey35', lwd = 0.25)
    
    # Print sample name at top left
    text(bp[1], ymax * 1.08, 
         labels = colnames(catalog)[ncol(catalog)], 
         xpd = NA,
         cex = cex, font = 2, adj = c(0, 0))

    if (catalog.type == "counts") {
      # Write the mutation counts on top of graph
      # count.cex <- 0.5 * (par("cex.main") - 1) + 1
      count.cex <- cex
      for (i in 1:6) {
        j <- 16 + 16 * (i - 1)
        k <- 1 + 16 * (i - 1)
        text(bp[j], ymax * 1.20, labels = sum(catalog[k:(16 * i), ]),
             adj = c(1, 1), xpd = NA, cex = count.cex)
      }
    }

    # Get locations for y axis annotations
    y.axis.values <- seq(0, ymax, ymax/4)
    if (catalog.type != "counts") {
      y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
    } else {
      y.axis.labels <- y.axis.values
    }
    
    if (grid) {
      segments(bp[1] - 0.5, seq(ymax/4, ymax, ymax/4), bp[num.classes] + 0.5,
               seq(ymax/4, ymax, ymax/4), col = 'grey35', lwd = 0.25)
      text(-0.5, y.axis.values, labels = y.axis.labels,
           las = 1, adj = 1, xpd = NA, cex = cex)
    }  else {
      Axis(side = 2, at = y.axis.values, las = 1, cex.axis = cex, labels = FALSE)
      text(-3.5, y.axis.values, labels = y.axis.labels, cex = cex,
           las = 1, adj = 1, xpd = NA)
    }
    
    if (xlabels) {
      # Draw the labels along bottom of x axis
      # cex.smaller <- cex * 0.8
      cex.xlabel <- cex
      xlabel.idx <- seq(1, 96, by = 4)
      label <- c("A", "C", "G", "T")

      # Draw the first line of x axis label
      text(bp[xlabel.idx], -ymax / 7, labels = label,
           cex = cex.xlabel, adj = 0.5, xpd = NA)

      x <- list(bp[xlabel.idx], bp[xlabel.idx + 1],
                bp[xlabel.idx + 2], bp[xlabel.idx + 3])
      y <- c(-ymax / 3.5, -ymax / 2.8, -ymax / 2.39, -ymax / 2.1)
      # Draw the remaining lines of x axis labels
      for (i in 1 : 4) {
        text(x[[i]], y[i], 
             labels = label[i], cex = cex.xlabel 
             , adj = 0.5, xpd = NA)
      }
      # Draw the text on the left plane
      text(1.5, -ymax / 7, labels = "preceded by 5'",
           pos = 2, xpd = NA, cex = cex.xlabel)
      text(1.5, -ymax / 3.5, labels = "followed by 3'",
           pos = 2, xpd = NA, cex = cex.xlabel)
    } else {
      every.fourth <- seq(from = 1, to = length(bp), by = 4)
      Axis(at = bp[every.fourth], side = 1, labels = FALSE, col = "grey35")
    }

    if (upper) {
      # Draw horizontal lines and names of major mutation class on top of graph
      x.left <- bp[seq(1, 81, 16)]
      x.right <- bp[seq(16, 96, 16)]
      rect(xleft = x.left, ybottom = ymax * 1.28, 
           xright = x.right, ytop = ymax * 1.3,
           col = class.col, border = NA, xpd = NA, adj = 0.5)
      text((x.left + x.right)/2, ymax * 1.38, 
           labels = maj.class.names, xpd = NA, cex = cex * 1.25)
    }

    return(list(plot.success = TRUE, plot.object = bp))
  }

#' @export
PlotCatalogToPdf.SBS96Catalog <-
  function(catalog, file, plot.SBS12, cex = 0.8,
           grid = TRUE, upper = TRUE, xlabels = TRUE,
           ylim = NULL) {
    old.par.tck.value <- par("tck")
    # Setting the width and length for A4 size plotting
    grDevices::cairo_pdf(file, width = 8.2677, 
                         height = 11.6929, onefile = TRUE)
    par(tck = old.par.tck.value)
    # opar <- par(no.readonly = TRUE)
    
    n <- ncol(catalog)
    opar <- par(mfrow = c(8, 1), mar = c(4, 5.5, 2, 1), oma = c(1, 1, 2, 1))
    on.exit(par(opar))

    for (i in 1 : n) {
      cat <- catalog[, i, drop = FALSE]
      PlotCatalog(cat, cex = cex, grid = grid, upper = upper, 
                  xlabels = xlabels, ylim = ylim)
    }
    
    grDevices::dev.off()
    return(list(plot.success = TRUE))
  }

#' @export
PlotCatalog.SBS192Catalog <- 
  function(catalog, plot.SBS12 = FALSE, cex = par("cex"),
           grid, upper, xlabels, ylim) {
  stopifnot(dim(catalog) == c(192, 1))

  if (plot.SBS12 == FALSE) {
    class.col  <- c("#03bcee",
                    "#010101",
                    "#e32926",
                    "#999999",
                    "#a1ce63",
                    "#ebc6c4")

    bg.class.col  <- c("#DCF8FF",
                       "#E9E9E9",
                       "#FFC7C7",
                       "#F7F7F7",
                       "#E5F9DF",
                       "#F9E7E7")

    strand.col <- c("#394398",
                    "#e83020")

    # Sort data in plotting order
    cat <- catalog[to.reorder.SBS.192.for.plotting, 1, drop = FALSE]

    num.classes <- length(cat)
    maj.class.names = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    cols <- rep(strand.col, num.classes / 2)

    if (attributes(cat)$catalog.type == "counts") {
      # Set a minimum value for ymax to make the plot more informative
      ymax <- 4 * ceiling(max(max(cat[, 1]) * 1.3, 10) / 4)

      # Barplot: side by side
      mat <- matrix(cat[, 1], nrow = 2, ncol = num.classes / 2)
      bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax),
                    axes = FALSE, lwd = 3, xaxs = "i", 
                    border = NA, col = cols, xpd = NA, ylab = "counts",
                    cex.lab = cex * par("cex.lab") * 1.25)
    } else if (attributes(cat)$catalog.type %in%
               c("counts.signature", "density.signature")) {
      # Determine the y axis label
      yaxislabel <- ifelse(attributes(cat)$catalog.type == "counts.signature",
                           "counts proportion", "density proportion")
      # Get ylim
      ymax <- ifelse(max(cat[, 1]) * 1.3 > 1, 1, max(cat[, 1]) * 1.3)

      # Barplot: side by side
      mat <- matrix(cat[, 1], nrow = 2, ncol = num.classes / 2)
      bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax),
                    axes = FALSE, lwd = 3, xaxs = "i",
                    border = NA, col = cols, xpd = NA, ylab = yaxislabel,
                    cex.lab = cex * par("cex.lab") * 1.25)
    } else if (attributes(cat)$catalog.type == "density") {
      # Get the rate of mutations per million trinucleotides
      rate <- cat[, 1] * 1000000

      # Get ylim
      ymax <- max(rate) * 1.3

      # Barplot: side by side
      mat <- matrix(rate, nrow = 2, ncol = num.classes / 2)
      bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax),
                    axes = FALSE, lwd = 3, xaxs = "i",
                    border = NA, col = cols, xpd = NA, ylab = "mut/million",
                    cex.lab = cex * par("cex.lab") * 1.25)
    }

    # Draw lines above each class:
    x.left <- bp[seq(0, 160, 32) + 1] - 0.5
    x.right <- bp[seq(32, 192, 32)] + 0.5
    rect(xleft = x.left, ymax * 1.01, xright = x.right, ymax * 1.03,
         col = class.col, border = NA, xpd = NA)

    # Draw mutation class labels at the top of the figure:
    text((x.left + x.right)/2, ymax * 1.07, labels = maj.class.names,
         cex = cex, xpd = NA)

    # Draw background color
    rect(xleft = x.left - 0.5, 0, xright = x.right + 0.5, ymax,
         col = bg.class.col, border = 'grey90', lwd = 1.5)

    # Plot again
    barplot(mat, beside = TRUE, ylim = c(0, ymax),
            axes = FALSE, ann = FALSE, lwd = 3,
            border = NA, col = cols, xpd = NA, add = TRUE)

    # Draw grid lines
    segments(bp[1] - 1, seq(0, ymax, ymax/4), bp[num.classes] + 1,
             seq(0, ymax, ymax/4), col = 'grey35', lwd = 0.25)

    # Draw y axis and write mutation counts on top of graph(if applicable)
    y.axis.values <- seq(0, ymax, ymax/4)
    if (attributes(cat)$catalog.type != "counts") {
      y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
    } else {
      y.axis.labels <- y.axis.values

      # Write the mutation counts on top of graph
      for (i in 1 : 6) {
        j <- 32 + 32 * (i - 1)
        k <- 1 + 32 * (i - 1)
        text(bp[j], ymax * 0.92, labels = sum(cat[k : (32 * i), 1]),
             adj = c(1, 1), xpd = NA, cex = cex)
      }
    }
    text(-0.5, y.axis.values, labels = y.axis.labels,
         las = 1, adj = 1, xpd = NA, cex = cex)

    # Draw the x axis labels
    cex.xlabel <- 0.625 * cex
    context.pos <- (bp[seq(1, 191, 2)] + bp[seq(2, 192, 2)]) / 2
    xlabel.1 <- c("A", "C", "G", "T")
    xlabel.2 <- rep(c("A", "C", "G", "T"), each = 4)
    text(context.pos, -ymax / 100, labels = rep(xlabel.1, 24), cex = cex.xlabel,
         srt = 90, adj = 1, xpd = NA)
    text(context.pos, -ymax / 18, labels = rep(c("C", "T"), each = 48),
         cex = cex.xlabel, srt = 90, adj = 1, xpd = NA)
    text(context.pos, -ymax / 10, labels = rep(xlabel.2, 6),
         cex = cex.xlabel, srt = 90, adj = 1, xpd = NA)

    # Write the name of the sample
    text(1.5, ymax * 7 / 8, labels = colnames(cat), adj = 0, cex = cex, font = 2)
  } else {
    strand.col <- c('#394398',
                    '#e83020')

    # Sort data in plotting order
    cat <- catalog[to.reorder.SBS.192.for.plotting, 1, drop = FALSE]

    num.classes <- 12
    maj.class.names <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    cols <- rep(strand.col, num.classes / 2)

    if (attributes(cat)$catalog.type == "counts") {
      # Get the counts for each major mutation class
      counts <- cat[, 1]
      counts.strand <- integer(12)
      for (i in 1 : 6){
        counts.strand[2 * i - 1] <-
          sum(counts[seq(32 * (i - 1) + 1, by = 2, length.out = 16)])
        counts.strand[2 * i] <-
          sum(counts[seq(32 * (i - 1) + 2, by = 2, length.out = 16)])
      }

      # Get ylim
      ymax <- 4 * ceiling(max(max(counts.strand) * 1.3, 10) / 4)

      # Barplot: side by side
      mat <- matrix(counts.strand, nrow = 2, ncol = num.classes / 2)
      bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 5.5),
                    width = 0.3, xaxs = "i", yaxs = "i",
                    axes = FALSE, ylab = "counts",
                    border = NA, col = cols, xpd = NA, 
                    cex.lab = cex * par("cex.lab") * 1.25)
      
      # Perform binomial test
      if (IsBinomialTestApplicable(cat)) {
        colnames(mat) <- maj.class.names
        rownames(mat) <- c("transcribed", "untranscribed")
        
        # Calculate the proportion of pyrimidines on transcribed strand
        # which can be used as the hypothesized probability of success
        # in binomial test
        counts <- CalBaseCountsFrom3MerAbundance(attributes(cat)$abundance)
        prop.C <- counts["C"] / sum(counts["C"] + counts["G"])
        prop.T <- counts["T"] / sum(counts["T"] + counts["A"])
        props <- c(rep(prop.C, 3), rep(prop.T, 3))
        names(props) <- maj.class.names
        
        p.values <- numeric(6)
        names(p.values) <- maj.class.names
        
        for (type in maj.class.names) {
          htest <- binom.test(x = mat[, type], p = props[type], 
                              alternative = "two.sided")
          p.values[type] <- htest$p.value
        }
        p.values <- p.adjust(p.values, method = "bonferroni")
        list0 <- list()
        list0[[colnames(cat)]] <- p.values
        
        # Draw asterisks on top of graph if p-value is significant
        for (type in maj.class.names) {
          p.value <- p.values[type]
          if (p.value < 0.05) {
            colnames(bp) <- maj.class.names
            # Get the x coordinates of the line segment to be drawn
            x1 <- bp[1, type]
            x2 <- bp[2, type]
            
            # Get the y coordinates of the line segment to be drawn
            y1 <- y2 <- max(mat[, type]) + max(mat) * 0.03
            
            # Draw the line segment
            segments(x1, y1, x2, y2)
            
            # Draw the asterisk on top of line segment
            x3 <- mean(c(x1, x2))
            y3 <- y1 + max(mat) * 0.025
            label <- AssignNumberOfAsterisks(p.value)
            text(x3, y3, label)
          }
        }
      }
    } else if (attributes(cat)$catalog.type %in%
               c("counts.signature", "density.signature")) {
      # Determine the y axis label
      yaxislabel <- ifelse(attributes(cat)$catalog.type == "counts.signature",
                           "counts proportion", "density proportion")

      # Get the proportion for each major mutation class
      prop <- cat[, 1]
      prop.strand <- integer(12)
      for (i in 1 : 6){
        prop.strand[2 * i - 1] <-
          sum(prop[seq(32 * (i - 1) + 1, by = 2, length.out = 16)])
        prop.strand[2 * i] <-
          sum(prop[seq(32 * (i - 1) + 2, by = 2, length.out = 16)])
      }

      # Get ylim
      ymax <- ifelse(max(prop.strand) * 1.3 > 1, 1, max(prop.strand) * 1.3)

      # Barplot: side by side
      mat <- matrix(prop.strand, nrow = 2, ncol = num.classes / 2)
      bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 5.5),
                    width = 0.3, xaxs = "i", yaxs = "i",
                    axes = FALSE, ylab = yaxislabel,
                    border = NA, col = cols, xpd = NA,
                    cex.lab = cex * par("cex.lab") * 1.25)
    } else if (attributes(cat)$catalog.type == "density") {
      # Get the rate of mutations per million trinucleotides
      rates <- cat[, 1] * 1000000

      rates.strand <- integer(12)
      for (i in 1 : 6){
        rates.strand[2 * i - 1] <-
          sum(rates[seq(32 * (i - 1) + 1, by = 2, length.out = 16)])
        rates.strand[2 * i] <-
          sum(rates[seq(32 * (i - 1) + 2, by = 2, length.out = 16)])
      }

      # Get ylim
      ymax <- max(rates.strand) * 1.3

      # Barplot: side by side
      mat <- matrix(rates.strand, nrow = 2, ncol = num.classes / 2)
      bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 5.5),
                    width = 0.3, xaxs = "i", yaxs = "i",
                    axes = FALSE, ylab = "mut/million",
                    border = NA, col = cols, xpd = NA,
                    cex.lab = cex * par("cex.lab") * 1.25)
    }

    # Draw y axis
    y.axis.values <- seq(0, ymax, length.out = 5)
    if (attributes(cat)$catalog.type != "counts") {
      y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
    } else {
      y.axis.labels <- round(y.axis.values, 0)
    }
    Axis(side = 2, at = y.axis.values, las = 1, labels = FALSE)
    text(-0.25, y.axis.values, labels = y.axis.labels,
         las = 1, adj = 1, xpd = NA, cex = cex)

    # Draw x axis
    xlabel <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    text(bp[seq(1, 12, by = 2)] + 0.16, -ymax / 28,
         labels = xlabel, xpd = NA, cex = cex)

    # Add legend
    legend(bp[6], ymax * 0.95, fill = strand.col, border = "white",
           xpd = NA, bty = "n",
           legend = c("Transcribed", "Untranscribed"), cex = cex)

    # Draw the sample name information on top of graph
    text(bp[5], ymax * 1.02, labels = colnames(catalog), xpd = NA,
         font = 2, cex = cex, adj = c(0, 0))
  }
    
    # Check whether it is possible to return the p-values from binomial test
    if (isTRUE(plot.SBS12) && IsBinomialTestApplicable(catalog)) {
      return(list(plot.success = TRUE, p.values = list0))
    } else {
      return(list(plot.success = TRUE))
    }
}

#' @export
PlotCatalogToPdf.SBS192Catalog <-
  function(catalog, file, plot.SBS12 = FALSE, cex = 0.8, 
           grid, upper, xlabels, ylim) {
  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  n <- ncol(catalog)
  
  ifelse(plot.SBS12, 
         par(mfrow = c(4, 3), mar = c(2, 5, 2, 1), oma = c(2, 2, 2, 2)),
         par(mfrow = c(8, 1), mar = c(2, 4, 2, 2), oma = c(3, 2, 1, 1)))
  
  p.values <- NULL
  for (i in 1:n) {
    cat <- catalog[, i, drop = FALSE]
    list <- PlotCatalog(cat, plot.SBS12 = plot.SBS12, cex = cex)
    p.values <- c(p.values, list$p.values)
  }
  grDevices::dev.off()
  ifelse(is.null(p.values), return(list(plot.success = TRUE)),
         return(list(plot.success = TRUE, p.values = p.values)))
}

#' @export
PlotCatalog.SBS1536Catalog <-
  function(catalog, plot.SBS12, cex, grid, upper, xlabels, ylim) {
  stopifnot(dim(catalog) == c(1536, 1))

  # Define the bases and their colors in plot
  bases <- c("A", "C", "G", "T")
  base.cols <- c("forestgreen", "dodgerblue2", "black", "red")

  # Define the theme color for plotting
  theme.col <- "forestgreen"
  colPal <- grDevices::colorRampPalette(c("white", theme.col))

  DrawImage <- function(counts, colors) {
    image(1 : 16, 1 : 16, matrix(counts, 16, 16),
          col = colors, axes = FALSE, ann = FALSE, bty = "n")

    # Draw the gridlines
    segments(rep(0.5, 5), c(0.5, 4.5, 8.5, 12.5, 16.5),
             rep(16.5, 5), c(0.5, 4.5, 8.5, 12.5, 16.5), xpd = NA)
    segments(c(0.5, 4.5, 8.5, 12.5, 16.5), rep(0.5, 5),
             c(0.5, 4.5, 8.5, 12.5, 16.5), rep(16.5, 5), xpd = NA)
  }

  DrawAxisY <- function() {
    text(0, 16 : 1, rep(bases, each = 4), col = rep(base.cols, each = 4),
         srt = 90, font = 2, xpd = NA)
    text(-1, 16 : 1, rep(bases, 4), col = rep(base.cols, 4),
         srt = 90, font = 2, xpd = NA)
    text(-2.5, 8.5, "Preceding bases", srt = 90, cex = 1.5, xpd = NA)
    segments(rep(-1.3, 5), c(0.5, 4.5, 8.5, 12.5, 16.5), rep(0.5, 5),
             c(0.5, 4.5, 8.5, 12.5, 16.5), xpd = NA)
  }

  DrawAxisX <- function() {
    text(1 : 16, 17, rep(bases, each = 4), col=rep(base.cols, each = 4),
         font = 2, xpd = NA)
    text(1 : 16, 18, rep(bases, 4), col = rep(base.cols, 4),
         font = 2, xpd = NA)
    segments(c(0.5, 4.5, 8.5, 12.5, 16.5), rep(16.5, 5),
             c(0.5, 4.5, 8.5, 12.5, 16.5), rep(18.3, 5), xpd = NA)
  }

  cat <- data.frame(catalog)
  cat$main.types <-
    paste0(substr(rownames(cat), 3, 3), ">", substr(rownames(cat), 6, 6))
  cat1 <-
    stats::aggregate(cat[, 1], by = list(main.types = cat$main.types), FUN = sum)

  if (attributes(catalog)$catalog.type == "counts") {
    # Get the total counts for the six main mutation types
    main.types.counts <- cat1[, 2]
    names(main.types.counts) <- cat1$main.types
  } else if (attributes(catalog)$catalog.type  %in%
             c("counts.signature", "density.signature")) {
    # Get the total proportion for the six main mutation types
    main.types.prop <- cat1[, 2]
    names(main.types.prop) <- cat1$main.types
  }

  # Sort the catalog matrix in plotting order
  rates <- as.data.frame(catalog)
  mut.type <- rownames(catalog)
  rates$mut.type <- mut.type
  rates$ref2alt <- paste0(substr(mut.type, 3, 3), substr(mut.type, 6, 6))
  rates$minus2bs <- substr(mut.type, 1, 1)
  rates$minus1bs <- substr(mut.type, 2, 2)
  rates$plus1bs <- substr(mut.type, 4, 4)
  rates$plus2bs <- substr(mut.type, 5, 5)
  rates <- dplyr::arrange(rates, ref2alt, dplyr::desc(minus1bs),
                          dplyr::desc(minus2bs), plus1bs, plus2bs)
  plot.order <- rates$mut.type
  rates <- as.matrix(rates[, 1])
  rownames(rates) <- plot.order

  main.mut.type <-
    paste(substr(plot.order, 3, 3), substr(plot.order, 6, 6), sep = '>')
  main.types <- unique(main.mut.type)
  n.types <- length(main.types)  # max 6 types

  # Plot one sample on one page
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mfrow = c(2, 3), oma = c(1, 6, 1, 4), pty = "s")

  for (i in 1 : n.types) {
    main.type <- main.types[i]
    sub.rates <- rates[main.mut.type == main.type, 1, drop = FALSE]
    col.ref <- colPal(256)

    # Draw the 6 plots on page one by one
    if (i == 1) {
      par(mar = c(0, 1, 7.5, 1))
      DrawImage(sub.rates, col.ref)
      DrawAxisY()
      DrawAxisX()
      text(8.5, 19, main.type, cex = 1.5, xpd = NA)
      if (attributes(catalog)$catalog.type == "counts") {
        text(11.5, 19, paste0("(N=", main.types.counts[main.type], ")"), xpd = NA)
      } else if (attributes(catalog)$catalog.type %in%
                 c("counts.signature", "density.signature")) {
        text(11.5, 19,
             paste0("(", round(100 * main.types.prop[main.type], 1), "%)"),
             xpd = NA)
      }
    }

    if (i == 2) {
      par(mar = c(0, 1, 7.5, 1))
      DrawImage(sub.rates, col.ref)
      DrawAxisX()
      text(8.5, 19, main.type, cex = 1.5, xpd = NA)
      if (attributes(catalog)$catalog.type == "counts") {
        text(11.5, 19, paste0("(N=", main.types.counts[main.type], ")"), xpd = NA)
      } else if (attributes(catalog)$catalog.type %in%
                 c("counts.signature", "density.signature")) {
        text(11.5, 19,
             paste0("(", round(100 * main.types.prop[main.type], 1), "%)"),
             xpd = NA)
      }
      text(8.5, 20.5, colnames(catalog), cex = 1.5, xpd = NA)
    }

    if (i == 3) {
      par(mar = c(0, 1, 7.5, 1))
      DrawImage(sub.rates, col.ref)
      DrawAxisX()
      text(8.5, 19, main.type, cex = 1.5, xpd = NA)
      if (attributes(catalog)$catalog.type == "counts") {
        text(11.5, 19, paste0("(N=", main.types.counts[main.type], ")"), xpd = NA)
      } else if (attributes(catalog)$catalog.type %in%
                 c("counts.signature", "density.signature")) {
        text(11.5, 19,
             paste0("(", round(100 * main.types.prop[main.type], 1), "%)"),
             xpd = NA)
        legend <- ifelse(attributes(catalog)$catalog.type == "counts.signature",
                         "counts proportion", "density proportion")
        text(10.3, 20.3, paste0("%: ", legend), xpd = NA)
      }
      text(17.5, 17, '1bp 3\'', xpd = NA, cex = 1)
      text(17.5, 18, '2bp 3\'', xpd = NA, cex = 1)
    }

    if (i == 4) {
      par(mar = c(6, 1, 0, 1))
      DrawImage(sub.rates, col.ref)
      DrawAxisY()
      text(8.5, 17, main.type, cex = 1.5, xpd = NA)
      if (attributes(catalog)$catalog.type  == "counts") {
        text(11.5, 17, paste0("(N=", main.types.counts[main.type], ")"), xpd = NA)
      } else if (attributes(catalog)$catalog.type %in%
                 c("counts.signature", "density.signature")) {
        text(11.5, 17,
             paste0("(", round(100 * main.types.prop[main.type], 1), "%)"),
             xpd = NA)
      }
      text(-1, -0.9, '1bp 5\'', xpd = NA, srt = 45, adj = 0, cex = 1)
      text(-2, -0.9, '2bp 5\'', xpd = NA, srt = 45, adj = 0, cex = 1)
    }

    if (i == 5) {
      par(mar = c(6, 1, 0, 1))
      DrawImage(sub.rates, col.ref)
      text(8.5, 17, main.type, cex = 1.5, xpd = NA)
      if (attributes(catalog)$catalog.type  == "counts") {
        text(11.5, 17, paste0("(N=", main.types.counts[main.type], ")"), xpd = NA)
      } else if (attributes(catalog)$catalog.type %in%
                 c("counts.signature", "density.signature")) {
        text(11.5, 17,
             paste0("(", round(100 * main.types.prop[main.type], 1), "%)"),
             xpd = NA)
      }
    }

    if (i == 6) {
      par(mar = c(6, 1, 0, 1))
      DrawImage(sub.rates, col.ref)
      text(8.5, 17, main.type, cex = 1.5, xpd = NA)
      if (attributes(catalog)$catalog.type  == "counts") {
        text(11.5, 17, paste0("(N=", main.types.counts[main.type], ")"), xpd = NA)
      } else if (attributes(catalog)$catalog.type %in%
                 c("counts.signature", "density.signature")) {
        text(11.5, 17,
             paste0("(", round(100 * main.types.prop[main.type], 1), "%)"),
             xpd = NA)
      }
    }
  }

  return(list(plot.success = TRUE))
}

#' @export
PlotCatalogToPdf.SBS1536Catalog <-
  function(catalog, file, plot.SBS12, cex, grid, upper, xlabels, ylim) {
  grDevices::cairo_pdf(file, width = 11.6929, height = 9.2677, onefile = TRUE)

  n <- ncol(catalog)

  for (i in 1 : n) {
    cat <- catalog[, i, drop = FALSE]
    PlotCatalog(cat)
  }

  grDevices::dev.off()
  return(list(plot.success = TRUE))
}

###############################################################################
# Plotting functions for DBS78, DBS144 and DBS136 catalog start here
###############################################################################

#' @export
PlotCatalog.DBS78Catalog <- function(catalog, plot.SBS12, cex,
                                     grid, upper, xlabels, ylim) {
  stopifnot(dim(catalog) == c(78, 1))
  stopifnot(rownames(catalog) == ICAMS::catalog.row.order$DBS78)

  dinuc.class.col <- RColorBrewer::brewer.pal(10, "Paired")
  cols <- rep(dinuc.class.col, c(9, 6, 9, 6, 9, 6, 6, 9, 9, 9))
  num.classes <- length(catalog)
  maj.class.names <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")

  if (attributes(catalog)$catalog.type == "density") {
    # Calculate rate of mutations per million nucleotides for the catalog
    rate <- catalog[, 1] * 1000000

    # Get ylim
    ymax <- ifelse(max(rate) * 1.3 > 1, 1, max(rate) * 1.3)

    # Barplot
    bp <- barplot(rate, ylim = c(0, ymax), xaxt = "n", yaxt = "n", xaxs = "i",
                  lwd = 3, space = 1.35, border = NA, col = cols,
                  xpd = NA, ylab = "mut/million", cex.lab = 0.8)
  } else if (attributes(catalog)$catalog.type == "counts") {
    # Set a minimum value for ymax to make the plot more informative
    ymax <- 4 * ceiling(max(max(catalog[, 1]) * 1.3, 10) / 4)

    # Barplot
    bp <- barplot(catalog[, 1], xaxt = "n", yaxt = "n", ylim = c(0, ymax),
                  xaxs = "i", lwd = 3, space = 1.35, border = NA, xaxs = "i",
                  col = cols, ylab = "counts", cex.lab = 0.8)

    # Write the mutation counts on top of graph
    for (i in 1 : 10) {
      j <- c(9, 15, 24, 30, 39, 45, 51, 60, 69, 78)
      name <- substr(rownames(catalog), 1, 2)
      text(bp[j[i] + 0.5], ymax * 0.92, xpd = NA, cex = 0.8,
           adj = c(1, 1), labels = sum(catalog[name == maj.class.names[i], ]))
    }
  } else if (attributes(catalog)$catalog.type %in%
             c("counts.signature", "density.signature")) {
    # Determine the y axis label
    yaxislabel <- ifelse(attributes(catalog)$catalog.type == "counts.signature",
                         "counts proportion", "density proportion")

    # Get ylim
    ymax <- ifelse(max(catalog[, 1]) * 1.3 > 1, 1, max(catalog[, 1]) * 1.3)

    # Barplot
    bp <- barplot(catalog[, 1], xaxt = "n", yaxt = "n", ylim = c(0, ymax),
                  lwd = 3, space = 1.35, border = NA, xaxs = "i",
                  col = cols, ylab = yaxislabel, cex.lab = 0.8)
  }

  # Draw box and grid lines
  rect(xleft = bp[1] - 1, 0, xright = bp[num.classes] + 1, ymax,
       border = "grey60", lwd = 0.5, xpd = NA)
  segments(bp[1] - 1, seq(0, ymax, ymax / 4), bp[num.classes] + 1,
           seq(0, ymax, ymax / 4), col = "grey60", lwd = 0.5, xpd = NA)

  # Draw lines above each class:
  x.left <- bp[c(0, 9, 15, 24, 30, 39, 45, 51, 60, 69) + 1] - 0.5
  x.right <- bp[c(9, 15, 24, 30, 39, 45, 51, 60, 69, 78)] + 0.5
  rect(xleft = x.left, ymax * 1.01, xright = x.right, ymax * 1.08,
       col = dinuc.class.col, border = NA, xpd = NA)

  # Draw mutation class labels at the top of the graph
  text((x.left + x.right) / 2, ymax * 1.123,
       labels = paste(maj.class.names, "NN", sep = ">"), cex = 0.7, xpd = NA)

  # Draw the sample name information on top of graph
  text(1.5, ymax * 7 / 8, labels = colnames(catalog), adj = 0, cex = 0.8, font = 2)

  # Draw y axis
  y.axis.values <- seq(0, ymax, ymax / 4)
  if (attributes(catalog)$catalog.type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  } else {
    y.axis.labels <- y.axis.values
  }
  text(0.35, y.axis.values, labels = y.axis.labels,
       las = 1, adj = 1, xpd = NA, cex = 0.75)

  # Draw x axis labels
  text(bp, -ymax / 80, labels = substr(rownames(catalog), 4, 4),
       cex = 0.5, srt = 90, adj = 1, xpd = NA)
  text(bp, -ymax / 15, labels = substr(rownames(catalog), 3, 3),
       cex = 0.5, srt = 90, adj = 1, xpd = NA)

  return(list(plot.success = TRUE))
}

#' @export
PlotCatalogToPdf.DBS78Catalog <-
  function(catalog, file, plot.SBS12, cex, grid, upper, xlabels, ylim) {
  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)

  # opar <- par(no.readonly = TRUE)
  n <- ncol(catalog)
  opar <- par(mfrow = c(8, 1), mar = c(2, 4, 2, 2), oma = c(3, 3, 2, 2))
  on.exit(par(opar))
  
  for (i in 1 : n) {
    cat <- catalog[, i, drop = FALSE]
    PlotCatalog(cat)
  }

  grDevices::dev.off()
  return(list(plot.success = TRUE))
}

#' @export
PlotCatalog.DBS144Catalog <- function(catalog, plot.SBS12, cex = par("cex"),
                                      grid, upper, xlabels, ylim) {
  stopifnot(dim(catalog) == c(144, 1))
  strand.col <- c('#394398',
                  '#e83020')

  # Sort data in plotting order
  cat <- catalog[to.reorder.DBS.144.for.plotting, 1, drop = FALSE]

  num.classes <- 20
  # maj.class.names <-
  #  c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  cols <- rep(strand.col, num.classes / 2)

  if (attributes(catalog)$catalog.type == "counts") {
    # Get the counts for each major mutation class
    counts <- cat[, 1]
    counts.strand <- integer(20)
    for (i in 1 : 10){
      idx <- c(0, 18, 24, 42, 48, 66, 72, 78, 96, 114, 132)
      counts.strand[2 * i - 1] <-
        sum(counts[seq(idx[i] + 1, idx[i + 1], by = 2)])
      counts.strand[2 * i] <-
        sum(counts[seq(idx[i] + 2, idx[i + 1], by = 2)])
    }
    
    # Set a minimum value for ymax to make the plot more informative
    ymax <- 4 * ceiling(max(max(counts.strand) * 1.3, 10) / 4)

    # Barplot: side by side
    mat <- matrix(counts.strand, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 9),
                  width = 0.3, xaxs = "i", yaxs = "i",
                  axes = FALSE, ylab = "counts",
                  border = NA, col = cols, xpd = NA,
                  cex.lab = cex * par("cex.lab"))
  } else if (attributes(catalog)$catalog.type %in%
             c("counts.signature", "density.signature")) {
    # Determine the y axis label
    yaxislabel <- ifelse(attributes(catalog)$catalog.type == "counts.signature",
                         "counts proportion", "density proportion")

    # Get the proportion for each major mutation class
    prop <- cat[, 1]
    prop.strand <- integer(20)
    for (i in 1 : 10){
      idx <- c(0, 18, 24, 42, 48, 66, 72, 78, 96, 114, 132)
      prop.strand[2 * i - 1] <-
        sum(prop[seq(idx[i] + 1, idx[i + 1], by = 2)])
      prop.strand[2 * i] <-
        sum(prop[seq(idx[i] + 2, idx[i + 1], by = 2)])
    }

    # Get ylim
    ymax <- ifelse(max(prop.strand) * 1.3 > 1, 1, max(prop.strand) * 1.3)

    # Barplot: side by side
    mat <- matrix(prop.strand, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 9),
                  width = 0.3, xaxs = "i", yaxs = "i",
                  axes = FALSE, ylab = yaxislabel,
                  border = NA, col = cols, xpd = NA,
                  cex.lab = cex * par("cex.lab"))
  } else if (attributes(catalog)$catalog.type == "density") {
    # Get the rate of mutations per million dinucleotides for each major mutation class
    rate <- cat[, 1] * 1000000
    rate.strand <- integer(20)
    for (i in 1 : 10){
      idx <- c(0, 18, 24, 42, 48, 66, 72, 78, 96, 114, 132)
      rate.strand[2 * i - 1] <-
        sum(rate[seq(idx[i] + 1, idx[i + 1], by = 2)])
      rate.strand[2 * i] <-
        sum(rate[seq(idx[i] + 2, idx[i + 1], by = 2)])
    }

    # Get ylim
    ymax <- max(rate.strand) * 1.3

    # Barplot: side by side
    mat <- matrix(rate.strand, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 9),
                  width = 0.3, xaxs = "i", yaxs = "i",
                  axes = FALSE, ylab = "mut/million",
                  border = NA, col = cols, xpd = NA,
                  cex.lab = cex * par("cex.lab"))
  }

  # Draw y axis
  y.axis.values <- seq(0, ymax, length.out = 5)
  if (attributes(catalog)$catalog.type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  } else {
    y.axis.labels <- y.axis.values
  }
  Axis(side = 2, at = y.axis.values, las = 1, labels = FALSE)
  text(-0.35, y.axis.values, labels = y.axis.labels,
       las = 1, adj = 1, xpd = NA, cex = cex)

  # Draw x axis
  xlabel <- c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  text(bp[seq(1, 20, by = 2)] + 0.16, -ymax / 28,
       labels = xlabel, xpd = NA, cex = cex)

  # Add legend
  legend(bp[10], ymax * 0.95, fill = strand.col, border = "white",
         xpd = NA, bty = "n",
         legend = c("Transcribed", "Untranscribed"), cex = cex)

  # Draw the sample name information on top of graph
  text(bp[8], ymax, labels = colnames(catalog), xpd = NA,
       font = 2, cex = cex, adj = c(0, 0))

  return(list(plot.success = TRUE))
}

#' @export
PlotCatalogToPdf.DBS144Catalog <-
  function(catalog, file, plot.SBS12, cex = 1, grid, upper, xlabels, ylim) {
    # Setting the width and length for A4 size plotting
    grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
    
    # opar <- par(no.readonly = TRUE)
    n <- ncol(catalog)
    opar <- par(mfrow = c(4, 3), mar = c(2, 5, 2, 1), oma = c(2, 2, 2, 2))
    on.exit(par(opar))
    
    for (i in 1 : n) {
      cat <- catalog[, i, drop = FALSE]
      PlotCatalog(cat, cex = cex)
    }

    grDevices::dev.off()
    return(list(plot.success = TRUE))
  }

#' @export
PlotCatalog.DBS136Catalog <- function(catalog, plot.SBS12, cex,
                                      grid, upper, xlabels, ylim) {
  stopifnot(dim(catalog) == c(136, 1))

  # Specify the layout of the plotting
  # invisible(
    layout(matrix(c(7, 8, 9, 10, 4, 5, 6, 11, 1, 2 , 3, 11), 3, 4,
                          byrow = TRUE))
  #  )

  # Define the bases and their colors in plot
  base <- c("A", "C", "G", "T")
  base.cols <- c("forestgreen", "dodgerblue2", "black", "red")

  ref.order <- c("AC", "AT", "GC", "CC", "CG", "CT", "TA", "TC", "TG", "TT")
  mut.type <- paste(ref.order, "NN", sep = ">")

  if (attributes(catalog)$catalog.type == "counts") {
    # Calculate the occurrences of each mutation type for plotting
    counts <- matrix(0, nrow = 160, ncol = 1)
    rownames(counts) <- order.for.DBS.136.plotting
    for (i in 1:160){
      if (order.for.DBS.136.plotting[i] %in% rownames(catalog)) {
        counts[i] <- catalog[order.for.DBS.136.plotting[i], ]
      } else {
        counts[i] <- NA
      }
    }

    # Calculate maximum count and total counts per mutation class
    df <- data.frame(stats::na.omit(counts))
    colnames(df) <- "counts"
    df$Ref <- substr(rownames(df), 2, 3)
    df1 <- stats::aggregate(df$counts, by = list(Ref = df$Ref), FUN = max)
    df2 <- stats::aggregate(df$counts, by = list(Ref = df$Ref), FUN = sum)
    max.count.per.class <- matrix(df1$x, 10, 1)
    counts.per.class <- matrix(df2$x, 10, 1)
    rownames(max.count.per.class) <- df1$Ref
    rownames(counts.per.class) <- df2$Ref
  } else if (attributes(catalog)$catalog.type == "density") {
    # Calculate tetranucleotide sequence contexts, normalized by tetranucleotide
    # occurrence in the genome
    rates <- matrix(0, nrow = 160, ncol = 1)
    rownames(rates) <- order.for.DBS.136.plotting
    for (i in 1:160){
      if (order.for.DBS.136.plotting[i] %in% rownames(catalog)) {
        rates[i] <- catalog[order.for.DBS.136.plotting[i], ]
      } else {
        rates[i] <- NA
      }
    }

    # Calculate maxima per mutation class(mut/million)
    df3 <- data.frame(stats::na.omit(rates))
    colnames(df3) <- "rates"
    df3$Ref <- substr(rownames(df3), 2, 3)
    df4 <- stats::aggregate(df3$rates, by = list(Ref = df3$Ref), FUN = max)
    max.rate.per.class <- matrix(round(df4$x * 1000000, 3), 10, 1)
    rownames(max.rate.per.class) <- df4$Ref
  } else {
    stop ('Plotting for DBS136 catlaog with type "', attributes(catalog)$catalog.type,
          '" is not implemented at the current stage.')
  }

  DrawImage <- function(mat) {
    value <- as.numeric(unlist(mat))
    maximum <- max(value[!is.na(value)])
    if (maximum == 0) {
      col.ref <- "white"
    } else {
      col.ref <- "forestgreen"
    }

    image(1:4, 1:4, mat,
          col = grDevices::colorRampPalette(c("white", col.ref))(16),
          axes = FALSE, ann = FALSE)

    # Make the background of the plot grey
    rect(0.5, 0.5, 4.5, 4.5 , col = "grey")

    # Plot the image again
    image(1:4, 1:4, mat,
          col = grDevices::colorRampPalette(c("white", col.ref))(16),
          axes = FALSE, ann = FALSE, add = TRUE)
  }

  opar <- par(mar = c(1, 2, 2, 0))
  on.exit(par(opar))
  for (i in 1:10){
    # par(mar = c(1, 2, 2, 0))

    if (attributes(catalog)$catalog.type == "density") {
      DrawImage(matrix(rates[(16 * (i - 1) + 1) : (16 * i)], 4, 4))
    } else if (attributes(catalog)$catalog.type == "counts") {
      DrawImage(matrix(counts[(16 * (i - 1) + 1) : (16 * i)], 4, 4))
    }

    # Draw the mutation type and number of occurrences on top of image
    text(2, 5.2, mut.type[i], font = 2, xpd = NA)
    if (attributes(catalog)$catalog.type == "counts") {
      text(3.2, 5.2, paste0("(", counts.per.class[ref.order[i], ], ")"), font = 2, xpd = NA)
    }

    # Draw a box surrounding the image
    segments(c(0.5, 0.5), c(0.5, 4.5), c(4.5, 4.5), c(0.5, 4.5), xpd = NA)
    segments(c(0.5, 4.5), c(0.5, 0.5), c(0.5, 4.5), c(4.5, 4.5), xpd = NA)

    # Draw the base information of the plot
    text(rep(0.2, 4), c(4, 3, 2, 1), base, col = base.cols, xpd = NA)
    text(seq(4), rep(4.8, 4), base, col = base.cols, xpd = NA)

    # Draw the sample name information of the sample
    if (i == 8) {
      mtext(colnames(catalog), at = 5, line =3)
    }
  }

  # Add in additional information
  plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
  text(x = 0.5, y = 0.9, "Maxima per class", cex = 1.6)
  ref <- c("TA", "TC", "TG", "TT", "CC", "CG", "CT", "AC", "AT", "GC")

  if (attributes(catalog)$catalog.type == "density") {
    text(x = 0.5, y = 0.8, "(mut/million)", cex = 1.2)
    maxima <- numeric(0)
    for (i in 1:10) {
      maxima[i] <- max.rate.per.class[ref[i], ]
      names(maxima)[i] <- ref[i]
    }
  } else if (attributes(catalog)$catalog.type == "counts") {
    text(x = 0.5, y = 0.8, "(counts)", cex = 1.2)
    maxima <- numeric(0)
    for (i in 1:10) {
      maxima[i] <- max.count.per.class[ref[i], ]
      names(maxima)[i] <- ref[i]
    }
  }

  text(rep(-0.1, 5), seq(0.7, 0.3, length.out = 5),
       paste(ref[1:5], maxima[1:5], sep = " = "), adj = 0, cex = 1.2, xpd = NA)
  text(rep(0.5, 5), seq(0.7, 0.3, length.out = 5),
       paste(ref[6:10], maxima[6:10], sep = " = "), adj = 0, cex = 1.2, xpd = NA)
  
  return(list(plot.success = TRUE))
}

#' @export
PlotCatalogToPdf.DBS136Catalog <-
  function(catalog, file, plot.SBS12, cex, grid, upper, xlabels, ylim) {
  stopifnot(nrow(catalog) == 136)
  n <- ncol(catalog)

  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(oma = c(1, 2, 1, 1))
  
  # Layout for plotting
  layout(matrix(c(12, 12, 12, 12,
                  7, 8, 9, 10, 7, 8, 9, 10,
                  7, 8, 9, 10, 7, 8, 9, 10,
                  4, 5, 6, 11,  4, 5, 6, 11,
                  4, 5, 6, 11, 4, 5, 6, 11,
                  1, 2, 3, 11, 1, 2, 3, 11,
                  1, 2, 3, 11,  1, 2, 3, 11,
                  24, 24, 24, 24,
                  19, 20, 21, 22, 19, 20, 21, 22,
                  19, 20, 21, 22, 19, 20, 21, 22,
                  16, 17, 18, 23, 16, 17, 18, 23,
                  16, 17, 18, 23, 16, 17, 18, 23,
                  13, 14, 15, 23,  13, 14, 15, 23,
                  13, 14, 15, 23, 13, 14, 15, 23),
                26, 4,byrow = TRUE))
  
  # Define the bases and their colors in plot
  base <- c("A", "C", "G", "T")
  base.cols <- c("forestgreen", "dodgerblue2", "black", "red")
  
  ref.order <- c("AC", "AT", "GC", "CC", "CG", "CT", "TA", "TC", "TG", "TT")
  mut.type <- paste(ref.order, "NN", sep = ">")
  
  for (i in 1:n) {
    cat <- catalog[, i, drop = FALSE]
    
    if (attributes(catalog)$catalog.type == "counts") {
      # Calculate the occurrences of each mutation type for plotting
      counts <- matrix(0, nrow = 160, ncol = 1)
      rownames(counts) <- order.for.DBS.136.plotting
      for (j in 1:160){
        if (order.for.DBS.136.plotting[j] %in% rownames(cat)) {
          counts[j] <- cat[order.for.DBS.136.plotting[j], ]
        } else {
          counts[j] <- NA
        }
      }
      
      # Calculate maximum count and total counts per mutation class
      df <- data.frame(stats::na.omit(counts))
      colnames(df) <- "counts"
      df$Ref <- substr(rownames(df), 2, 3)
      df1 <- stats::aggregate(df$counts, by = list(Ref = df$Ref), FUN = max)
      df2 <- stats::aggregate(df$counts, by = list(Ref = df$Ref), FUN = sum)
      max.count.per.class <- matrix(df1$x, 10, 1)
      counts.per.class <- matrix(df2$x, 10, 1)
      rownames(max.count.per.class) <- df1$Ref
      rownames(counts.per.class) <- df2$Ref
    }
    
    if (attributes(catalog)$catalog.type == "density") {
      # Calculate tetranucleotide sequence contexts, normalized by tetranucleotide
      # occurrence in the genome
      rates <- matrix(0, nrow = 160, ncol = 1)
      rownames(rates) <- order.for.DBS.136.plotting
      for (j in 1:160){
        if (order.for.DBS.136.plotting[j] %in% rownames(cat)) {
          rates[j] <- cat[order.for.DBS.136.plotting[j], ]
        } else {
          rates[j] <- NA
        }
      }
      
      # Calculate maxima per mutation class(mut/million)
      df3 <- data.frame(stats::na.omit(rates))
      colnames(df3) <- "rates"
      df3$Ref <- substr(rownames(df3), 2, 3)
      df4 <- stats::aggregate(df3$rates, by = list(Ref = df3$Ref), FUN = max)
      max.rate.per.class <- matrix(round(df4$x * 1000000, 3), 10, 1)
      rownames(max.rate.per.class) <- df4$Ref
    }

    DrawImage <- function(mat) {
      value <- as.numeric(unlist(mat))
      maximum <- max(value[!is.na(value)])
      if (maximum == 0) {
        col.ref <- "white"
      } else {
        col.ref <- "forestgreen"
      }

      image(1:4, 1:4, mat,
            col = grDevices::colorRampPalette(c("white", col.ref))(16),
            axes = FALSE, ann = FALSE)

      # Make the background of the plot grey
      rect(0.5, 0.5, 4.5, 4.5 , col = "grey")

      # Plot the image again
      image(1:4, 1:4, mat,
            col = grDevices::colorRampPalette(c("white", col.ref))(16),
            axes = FALSE, ann = FALSE, add = TRUE)
    }

    for (j in 1:10) {
      par(mar = c(1, 0, 2, 0), pty = "s")
      if (attributes(catalog)$catalog.type == "density") {
        DrawImage(matrix(rates[(16 * (j - 1) + 1) : (16 * j)], 4, 4))
      } else if (attributes(catalog)$catalog.type == "counts") {
        DrawImage(matrix(counts[(16 * (j - 1) + 1) : (16 * j)], 4, 4))
      } else {
        stop('Please specify the correct type: "density" or "counts"')
      }

      # Draw the mutation type and number of occurrences on top of image
      text(2.3, 5.2, mut.type[j], font = 2, xpd = NA)
      if (attributes(catalog)$catalog.type == "counts") {
        text(3.5, 5.2, paste0("(", counts.per.class[ref.order[j], ], ")"),
             font = 2, xpd = NA)

      }

      # Draw a box surrounding the image
      segments(c(0.5, 0.5), c(0.5, 4.5), c(4.5, 4.5), c(0.5, 4.5), xpd = NA)
      segments(c(0.5, 4.5), c(0.5, 0.5), c(0.5, 4.5), c(4.5, 4.5), xpd = NA)

      # Draw the base information of the plot
      text(rep(0.2, 4), c(4, 3, 2, 1), base, col = base.cols, xpd = NA)
      text(seq(4), rep(4.8, 4), base, col = base.cols, xpd = NA)
    }

    # Add in additional information
    plot.new()
    text(x = 0.5, y = 1.2, "Maxima per class", cex = 1.6, xpd = NA)
    ref <- c("TA", "TC", "TG", "TT", "CC", "CG", "CT", "AC", "AT", "GC")

    if (attributes(catalog)$catalog.type == "density") {
      text(x = 0.5, y = 1.07, "(mut/million)", cex = 1.2, xpd = NA)
      maxima <- numeric(0)
      for (j in 1:10) {
        maxima[j] <- max.rate.per.class[ref[j], ]
        names(maxima)[j] <- ref[j]
      }
    } else if (attributes(catalog)$catalog.type == "counts") {
      text(x = 0.5, y = 1.07, "(counts)", cex = 1.2, xpd = NA)
      maxima <- numeric(0)
      for (j in 1:10) {
        maxima[j] <- max.count.per.class[ref[j], ]
        names(maxima)[j] <- ref[j]
      }
    } else {
      stop('Please specify the correct type: "density" or "counts"')
    }
    text(rep(0, 5), seq(0.9, 0.1, length.out = 5),
         paste(ref[1:5], maxima[1:5], sep = " = "), adj = 0, cex = 1.2, xpd = NA)
    text(rep(0.6, 5), seq(0.9, 0.1, length.out = 5),
         paste(ref[6:10], maxima[6:10], sep = " = "), adj = 0, cex = 1.2, xpd = NA)

    # Draw the sample name information of the sample
    par(mar = c(0, 0, 0, 0))
    plot.new()
    text(0.7, 0.5, colnames(catalog)[i], cex = 1.5, xpd = NA)
  }

  grDevices::dev.off()
  return(list(plot.success = TRUE))
}

###############################################################################
# Plotting functions for ID(insertion and deletion) catalog start here
###############################################################################

#' @export
PlotCatalog.IndelCatalog <- function(catalog, plot.SBS12, cex,
                                     grid, upper, xlabels, ylim){
  stopifnot(dim(catalog) == c(83, 1))

  indel.class.col <- c("#fdbe6f",
                       "#ff8001",
                       "#b0dd8b",
                       "#36a12e",
                       "#fdcab5",
                       "#fc8a6a",
                       "#f14432",
                       "#bc141a",
                       "#d0e1f2",
                       "#94c4df",
                       "#4a98c9",
                       "#1764ab",
                       "#e2e2ef",
                       "#b6b6d8",
                       "#8683bd",
                       "#61409b")

  num.classes <- length(catalog)
  cols <- rep(indel.class.col,
              c(6, 6, 6, 6,
                6, 6, 6, 6,
                6, 6, 6, 6,
                1, 2, 3, 5))

  if (attributes(catalog)$catalog.type == "counts") {
    # Set a minimum value for ymax to make the plot more informative
    ymax <- 4 * ceiling(max(max(catalog[, 1]) * 1.3, 10) / 4)

    # Barplot
    bp <- barplot(catalog[, 1], ylim = c(0, ymax), axes = FALSE, xaxt = "n",
                  lwd = 3, space = 1.35, border = NA, col = cols, xpd = NA,
                  xaxs = "i", yaxt = "n")

    # Calculate and draw the total counts for each major type
    counts <- integer(16)
    for (i in 1:16) {
      idx <- c(6 * 1:12, 73, 75, 78, 83)
      idx2 <- c(8.9, 23, 37.1, 51.2, 65.3, 79.4, 93.5,
                107.6, 121.7, 135.8, 149.9, 164,
                172.2, 175.5, 182, 191)
      if (i == 1) {
        counts[i] <- sum(catalog[1:idx[1], 1])
      } else {
        counts[i] <- sum(catalog[(idx[i - 1] + 1):idx[i], 1])
      }
      text(idx2[i], ymax * 0.6, labels = counts[i],
           cex = 0.68, adj = 1, xpd = NA)
    }

  } else if (attributes(catalog)$catalog.type == "counts.signature") {
    # Get ylim
    ymax <- ifelse(max(catalog[, 1]) * 1.3 > 1, 1, max(catalog[, 1]) * 1.3)

    # Barplot
    bp <- barplot(catalog[, 1], ylim = c(0, ymax), axes = FALSE, xaxt = "n",
                  lwd = 3, space = 1.35, border = NA, col = cols, xpd = NA,
                  xaxs = "i", yaxt = "n")
  }

  # Draw box and grid lines
  rect(xleft = bp[1] - 1.5, 0, xright = bp[num.classes] + 1, ymax, col = NA,
       border = "grey60", lwd = 0.5, xpd = NA)
  segments(bp[1] - 1.5, seq(0, ymax, ymax / 4), bp[num.classes] + 1,
           seq(0, ymax, ymax / 4), col = "grey60", lwd = 0.5, xpd = NA)

  # Draw mutation class labels and lines above each class
  maj.class.names <- c("1bp deletion", "1bp insertion",
                       ">1bp deletions at repeats\n(Deletion length)",
                       ">1bp insertions at repeats\n(Insertion length)",
                       "Deletions with microhomology\n(Deletion length)")
  x.left <- bp[c(seq(0, 66, 6), 72, 73, 75, 78) + 1] - 0.5
  x.right <- bp[c(seq(6, 72, 6), 73, 75, 78, 83)] + 0.5
  class.pos <- c((x.left[seq(1, 4, 2)] + x.right[seq(2, 5, 2)]) / 2,
                 (x.left[c(6, 10)] + x.right[c(8, 12)] - 12) / 2,
                 (x.left[13] + x.right[length(x.left)]) / 2)
  category.lab <- c(rep(c("C", "T"), 2), rep(c("2", "3", "4", "5+"), 3))
  category.col <- c(rep(c("black", "white"), 2),
                    rep(c("black", "black", "black", "white"), 3))

  # Draw lines above each class
  rect(xleft = x.left, ymax * 1.02, xright = x.right, ymax * 1.11,
       col = indel.class.col, border = NA, xpd = NA)
  text((x.left + x.right) / 2, ymax * 1.06, labels = category.lab,
       cex = 0.65, col = category.col, xpd = NA)

  # Draw mutation class labels at the top of the figure
  text(class.pos, ymax * 1.27, labels = maj.class.names, cex = 0.75, xpd = NA)

  # Draw the sample name information of the sample
  text(1.5, ymax * 7 / 8, labels = colnames(catalog),
       adj = 0, cex = 0.8, font = 2)

  # Draw y axis
  y.axis.values <- seq(0, ymax, ymax / 4)
  if (attributes(catalog)$catalog.type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
    text(-9, ymax / 2, labels = "counts proportion",
         srt = 90, xpd = NA, cex = 0.8)
  } else {
    y.axis.labels <- y.axis.values
    text(-9, ymax / 2, labels = "counts",
         srt = 90, xpd = NA, cex = 0.8)
  }
  text(0, y.axis.values, labels = y.axis.labels,
       las = 1, adj = 1, xpd = NA, cex = 0.75)


  # Draw x axis labels
  mut.type <- c(rep(c("1", "2", "3", "4", "5", "6+"), 2),
                rep(c("0", "1", "2", "3", "4", "5+"), 2),
                rep(c("1", "2", "3", "4", "5", "6+"), 4),
                rep(c("0", "1", "2", "3", "4", "5+"), 4),
                "1", "1", "2", "1", "2", "3", "1", "2", "3", "4", "5+")
  bottom.pos <- c((x.left[1] + x.right[2]) / 2, (x.left[3] + x.right[4]) / 2,
                  class.pos[3 : length(class.pos)])
  bottom.lab <- c("Homopolymer length", "Homopolymer length",
                  "Number of repeat units", "Number of repeat units",
                  "Microhomology length")
  rect(xleft = x.left, -ymax * 0.09, xright = x.right, -ymax * 0.01,
       col = indel.class.col, border = NA, xpd = NA)
  text(bp, -ymax * 0.15, labels = mut.type, cex = 0.65, xpd = NA)
  text(bottom.pos, -ymax * 0.27, labels = bottom.lab, cex = 0.75, xpd = NA)

  return(list(plot.success = TRUE))
}

#' @export
PlotCatalogToPdf.IndelCatalog <-
  function(catalog, file, plot.SBS12, cex, grid, upper, xlabels, ylim) {
  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  
  n <- ncol(catalog)
  opar <- par(mfrow = c(8, 1), mar = c(3, 4, 2.5, 2), oma = c(3, 3, 2, 2))
  on.exit(par(opar))
  
  for (i in 1 : n) {
    cat <- catalog[, i, drop = FALSE]
    PlotCatalog(cat)
  }
  grDevices::dev.off()
  return(list(plot.success = TRUE))
}


#' Plot position probability matrix (PPM) for *one* sample from a Variant Call Format
#' (VCF) file.
#'
#' @param pwm A position probability matrix (PPM) for *one* sample.
#'
#' @param title The main title of the plot.
#'
#' @return \code{invisible(TRUE)}
#'
#' @keywords internal
PlotPPM <- function(ppm, title) {
  ppm <- t(ppm)
  x <- seq(1:nrow(ppm))
  plot(x, y = ppm[, "A"], xaxt = "n", xlab = "",
       ylab = "proportion", ylim = c(0, 1), pch = 20, type = "b",
       lwd = 2, col = "darkgreen")
  lines(x, ppm[, "C"], col = "blue", type = "b", pch = 20, lwd = 2)
  lines(x, ppm[, "G"], col = "black", type = "b", pch = 20, lwd = 2)
  lines(x, ppm[, "T"], col = "red", type = "b", pch = 20, lwd = 2)
  abline(h = 0.25, col = "grey50")
  text(x, y = -0.1, labels = rownames(ppm), adj = 1, srt = 90, xpd = NA)
  legend("topright", pch = 16, ncol = 4, bty = "n",
         legend = c("A", "C", "G", "T"),
         col = c("darkgreen", "blue", "black", "red"))
  mtext(text = title, side = 3, line = 0.5)
  invisible(TRUE)
}


#' Plot position probability matrices (PPM) to a PDF file
#'
#' @param list.of.pwm A list of position probability matrices (PPM)
#'
#' @param file The name of the PDF file to be produced.
#'
#' @param titles A vector of titles on top of each PPM plot.
#'
#' @return \code{invisible(TRUE)}
#'
#' @keywords internal
PlotPPMToPdf <- function(list.of.ppm, file, titles = names(list.of.ppm)) {
  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  
  n <- length(list.of.ppm)
  opar <- par(mfrow = c(4, 2), mar = c(4, 5.5, 2, 1), oma = c(1, 1, 2, 1))
  on.exit(par(opar))
  
  for (i in 1:n) {
    PlotPPM(list.of.ppm[[i]], title = titles[i])
  }
  
  grDevices::dev.off()
  invisible(TRUE)
}

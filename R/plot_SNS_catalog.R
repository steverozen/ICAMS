#' Plot Catalog Functions
#'
#' Plot the mutation catalog of one sample
#'
#' \code{PlotCat96} Plot the SNS 96 mutation catalog of one sample.
#'
#' \code{PlotCat192} Plot the SNS 192 mutation catalog of one sample.
#'
#' \code{PlotCat192Strand} Plot the transcription strand bias graph of 6 SNS
#' mutation types ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") in one sample.
#'
#' \code{PlotCat1536} Plot the pentanucleotide sequence contexts for one sample,
#' normalized by pentanucleotide occurrence in the genome. The mutation types
#' are in six-letters like CATTAT, first 2-letters CA refers to (-2, -1)
#' position, third letter T refers to the base which has mutation, next second
#' 2-letters TA refers to (+1, +2) position, last letter T refers to the base
#' after mutation.
#'
#' \code{PlotCatDNS78} Plot the DNS 78 mutation catalog of one sample.
#'
#' \code{PlotCatDNS144} Plot the transcription strand bias graph of 10 major DNS
#' mutation types ("AC>NN", "AT>NN", "CC>NN", "CG>NN", "CT>NN", "GC>NN",
#' "TA>NN", "TC>NN", "TG>NN", "TT>NN") in one sample.
#'
#' \code{PlotQUAD136} Plot the tetranucleotide sequence contexts of 10 major DNS
#' mutation types ("AC>NN", "AT>NN", "CC>NN", "CG>NN", "CT>NN", "GC>NN",
#' "TA>NN", "TC>NN", "TG>NN", "TT>NN") for one sample.
#'
#' \code{PlotCatID} Plot the insertion and deletion catalog of one sample.
#' (Please take note that the deletions Repeat Size ranges from 0 to 5+ in the
#' catalog, but for plotting and end user documentation it ranges from 1 to 6+.)
#' @param catalog A matrix whose rownames indicate the mutation types
#'   while its columns contain the counts of each mutation type.
#' @param id The ID information of the sample which has mutations.
#' @param type A value indicating the type of graph. If type = "counts", the
#'   graph will plot the occurrences of the mutation types in the sample. If
#'   type = "signature", the graph will plot mutation signatures of the sample.
#'   If type = "density", the graph will plot the rates of mutations per million
#'   nucleotides for each mutation type. (Please take note there is no "density"
#'   type for PlotCatID function and the option of type = "density" is not
#'   implemented for function PlotCat192, PlotCat192Strand and PlotCatDNS144 at
#'   the current stage.)
#' @param cex A numerical value giving the amount by which mutation class labels,
#'   mutation counts(if there exists), y axis and its labels, x axis labels and
#'   its annotations(if there exists) sample name and legend(if there exists)
#'   should be magnified relative to the default.
#' @param grid A logical value indicating whether to draw the grid lines in the
#'   graph.
#' @param upper A logical value indicating whether to draw horizontal lines and
#'   names of major mutation class on top of graph.
#' @param xlabels A logical value indicating whether to draw x axis labels.
#' @param abundance A matrix containing nucleotide abundance information and
#'   strand information(if there exists), to be used only when type = "density".
#' @return invisible(TRUE)
#' @name PlotCatalog
NULL


#' Catalog to Pdf Functions
#'
#' Plot the mutation catalog of different samples to a PDF file
#'
#' \code{Cat96ToPdf} Plot the SNS 96 mutation catalog of different samples
#' to a PDF file.
#'
#' \code{Cat192ToPdf} Plot the SNS 192 mutation catalog of different samples
#' to a PDF file.
#'
#' \code{Cat192StrandToPdf} Plot the transcription strand bias graph of
#' 6 SNS mutation types ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
#' of different samples to a PDF file.
#'
#' \code{Cat1536ToPdf} Plot the 1536 mutation catalog of >= 1 samples to a PDF
#' file. The mutation types are in six-letters like CATTAT, first 2-letters CA
#' refers to (-2, -1) position, third letter T refers to the base which has
#' mutation, next second 2-letters TA refers to (+1, +2) position, last letter T
#' refers to the base after mutation.
#'
#' \code{CatDNS78ToPdf} Plot the DNS 78 mutation catalog of different samples
#' to a PDF file.
#'
#' \code{CatDNS144ToPdf} Plot the transcription strand bias graph of
#' 10 major DNS mutation types ("AC>NN", "AT>NN", "CC>NN", "CG>NN", "CT>NN",
#' "GC>NN", "TA>NN", "TC>NN", "TG>NN", "TT>NN") of different samples
#' to a PDF file.
#'
#' \code{CatQUAD136ToPdf} Plot the tetranucleotide sequence contexts of 10 major
#' DNS mutation types ("AC>NN", "AT>NN", "CC>NN", "CG>NN", "CT>NN", "GC>NN",
#' "TA>NN", "TC>NN", "TG>NN", "TT>NN") of different samples to a PDF file.
#'
#' \code{CatIDToPdf} Plot the insertion and deletion catalog of different
#' samples to a PDF file. (Please take note that the deletions Repeat Size
#' ranges from 0 to 5+ in the catalog, but for plotting and end user
#' documentation it ranges from 1 to 6+.)
#' @param catalog A matrix whose rownames indicate the mutation types
#'   while its columns contain the counts of each mutation type from
#'   different samples.
#' @param name The name of the PDF file to be produced.
#' @param id A vector containing the ID information of different samples.
#' @param type A vector of values indicating the type of plot for each sample.
#'   If type = "counts", the graph will plot the occurrences of the mutation
#'   types in the sample. If type = "signature", the graph will plot mutation
#'   signatures of the sample. If type = "density", the graph will plot the
#'   rates of mutations per million nucleotides for each mutation type. (Please
#'   take note there is no "density" type for CatIDtoPdf function and the option
#'   of type = "density" is not implemented for function Cat192ToPdf,
#'   Cat192StrandToPdf and CatDNS144ToPdf at the current stage.)
#' @param cex A numerical value giving the amount by which mutation class labels,
#'   y axis labels, sample name and legend(if there exists) should be magnified
#'   relative to the default.
#' @param grid A logical value indicating whether to draw the grid lines in the
#'   graph.
#' @param upper A logical value indicating whether to draw horizontal lines and
#'   names of major mutation class on top of graph.
#' @param xlabels A logical value indicating whether to draw x axis labels.
#' @param abundance A matrix containing nucleotide abundance information, to
#'   be used only when type = "density".
#' @return invisible(TRUE)
#' @name CatalogToPdf
NULL

#' @rdname PlotCatalog
#' @import graphics
#' @export
PlotCat96 <-
  function(catalog, id = colnames(catalog), type = "density", cex = 0.8, grid = TRUE,
           upper = TRUE, xlabels = TRUE, abundance = NULL) {
    stopifnot(dim(catalog) == c(96, 1))
    stopifnot(rownames(catalog) == .catalog.row.order96)

    class.col <- c("#0000ff",  # dark blue
                   "#000000",  # black
                   "#ff4040",  # red
                   "#838383",  # grey
                   "#40ff40",  # green
                   "#ff667f")  # pink

    cols <- rep(class.col, each = 16)
    maj.class.names <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    num.classes <- length(catalog)

    if (type == "density") {
      # Calculate rate of mutations per million trinucleotides for the catalog
      rate <- double(96)
      for (i in 1 : 96) {
        rate[i] <-
          catalog[i] * 1000000 / abundance[substr(rownames(catalog)[i], 1, 3), ]
      }

      # Get ylim
      ymax <- max(rate)

      # Barplot
      bp <- barplot(rate, xaxt = "n", yaxt = "n", xaxs = "i",
                    xlim = c(-1, 230), lwd = 3, space = 1.35, border = NA,
                    col = cols, ylab = "mut/million")

      # Write the mutation counts on top of graph
      for (i in 1 : 6) {
        j <- 13 + 16 * (i - 1)
        k <- 1 + 16 * (i - 1)
        text(bp[j], ymax * 1.15, labels = sum(catalog[k : (16 * i), ]),
             xpd = NA, cex = cex)
      }
    } else if (type == "counts") {
      # Get ylim
      ymax <- max(catalog[, 1])

      # Barplot
      bp <- barplot(catalog[, 1], xaxt = "n", yaxt = "n", xlim = c(-1, 230),
                    xaxs = "i", lwd = 3, space = 1.35, border = NA,
                    col = cols, ylab = "counts")

      # Write the mutation counts on top of graph
      for (i in 1 : 6) {
        j <- 13 + 16 * (i - 1)
        k <- 1 + 16 * (i - 1)
        text(bp[j], ymax * 1.15, labels = sum(catalog[k : (16 * i), ]),
             xpd = NA, cex = cex)
      }
    } else if (type == "signature") {
      # Calculate mutation signatures of the input catalog
      sig <- catalog / sum(catalog)

      # Get ylim
      ymax <- max(sig)

      # Barplot
      bp <- barplot(sig[, 1], xaxt = "n", yaxt = 'n', xaxs = "i", xlim = c(-1, 230),
                    lwd = 3, space = 1.35, border = NA,
                    col = cols, ylab = "proportion")
    } else {
      stop('Please specify the correct type: "density", "counts" or "signature"')
    }

    # Draw grid lines?
    if (grid) {
      segments(bp[1] - 1, seq(ymax/4, ymax, ymax/4), bp[num.classes] + 1,
               seq(ymax/4, ymax, ymax/4), col = 'grey35', lwd = 0.25)
    }

    # Draw the x axis
    Axis(side = 1, at = c(bp[seq(1, 93, 4)] - 0.5, bp[96] + 0.5),
         labels = FALSE, lwd.tick = 0, lwd = 0.5)

    # Draw y axis
    y.axis.values <- seq(0, ymax, ymax/4)
    if (type != "counts") {
      y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
    } else {
      y.axis.labels <- round(y.axis.values, 0)
    }
    if (grid) {
      text(-0.5, y.axis.values, labels = y.axis.labels,
           las = 1, adj = 1, xpd = NA, cex = cex)
    } else {
      Axis(side = 2, at = y.axis.values, las = 1, cex.axis = cex, labels = FALSE)
      text(-3.5, y.axis.values, labels = y.axis.labels, cex = cex,
           las = 1, adj = 1, xpd = NA)
    }

    # Draw the ID information on top of graph
    text(bp[2], ymax * 1.08, labels = id, xpd = NA, font = 2, adj = c(0, 0))

    # Draw the labels along x axis?
    if (xlabels) {
      xlabel.idx <- seq(1, 96, by = 4)
      label <- c("A", "C", "G", "T")

      # Draw the first line of x axis label
      text(bp[xlabel.idx], -ymax / 7, labels = label,
           cex = cex, adj = 0.5, xpd = NA)

      x <- list(bp[xlabel.idx], bp[xlabel.idx + 1],
                bp[xlabel.idx + 2], bp[xlabel.idx + 3])
      y <- c(-ymax / 3.5, -ymax / 2.8, -ymax / 2.5, -ymax / 2.1)
      # Draw the remaining lines of x axis labels
      for (i in 1 : 4) {
        text(x[[i]], y[i], labels = label[i], cex = cex, adj = 0.5, xpd = NA)
      }
    }

    # Draw the text on the left plane
    text(1.5, -ymax / 7, labels = "preceded by 5'",
         pos = 2, xpd = NA, cex = cex)
    text(1.5, -ymax / 3.5, labels = "followed by 3'",
         pos = 2, xpd = NA, cex = cex)

    # Draw horizontal lines and names of major mutation class on top of graph?
    if (upper) {
      x.left <- bp[seq(1, 81, 16)]
      x.right <- bp[seq(16, 96, 16)]
      rect(xleft = x.left, ymax * 1.28, xright = x.right, ymax * 1.3,
           col = class.col, border = NA, xpd = NA, adj = 0.5)
      text((x.left + x.right)/2, ymax * 1.38, labels = maj.class.names, xpd = NA)
    }

    invisible(TRUE)
  }

#' @rdname CatalogToPdf
#' @export
Cat96ToPdf <-
  function(catalog, name, id = colnames(catalog), type = "density",
           grid = FALSE, upper = TRUE, xlabels = TRUE, abundance = NULL) {
    # Setting the width and length for A4 size plotting
    grDevices::cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE)

    n <- ncol(catalog)
    graphics::par(mfrow = c(8, 1), mar = c(4, 5.5, 2, 1), oma = c(1, 1, 2, 1))

    # Do recycling of the function parameters if a vector
    # with length more than one is not specified by the user.
    if (n > 1 && length(type) == 1) {
      type <- rep(type, n)
    }

    if (n > 1 && length(grid) == 1) {
      grid <- rep(grid, n)
    }

    if (n > 1 && length(upper) == 1) {
      upper <- rep(upper, n)
    }

    if (n > 1 && length(xlabels) == 1) {
      xlabels <- rep(xlabels, n)
    }

    for (i in 1 : n) {
      PlotCat96(catalog[, i, drop = FALSE],
                id = id[i],
                type = type[i],
                grid = grid[i],
                upper = upper[i],
                xlabels = xlabels[i],
                abundance = abundance)
    }
    invisible(grDevices::dev.off())
    invisible(TRUE)
  }

#' @rdname PlotCatalog
#' @import graphics
#' @export
PlotCat192 <- function(catalog, id = colnames(catalog), type = "counts",
                       cex = 0.8, abundance = NULL) {
  stopifnot(dim(catalog) == c(192, 1))

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
  counts <- catalog[.to.reorder.192.for.plotting, ]

  num.classes <- length(counts)
  maj.class.names = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  cols <- rep(strand.col, num.classes / 2)

  if (type == "counts") {
    # Get ylim
    ymax <- max(counts) * 1.3

    # Barplot: side by side
    mat <- matrix(counts, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax),
                  axes = FALSE, ann = FALSE, lwd = 3, xaxs = "i",
                  border = NA, col = cols, xpd = NA, ylab = "counts")
  } else if (type == "signature") {
    # Calculate mutation signatures of the input catalog
    sig <- counts / sum(counts)

    # Get ylim
    ymax <- ifelse(max(sig) * 1.3 > 1, 1, max(sig) * 1.3)

    # Barplot: side by side
    mat <- matrix(sig, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax),
                  axes = FALSE, ann = FALSE, lwd = 3, xaxs = "i",
                  border = NA, col = cols, xpd = NA, ylab = "proportion")
  } else if (type == "density") {
    stop('type = "density" not implemented')
  } else {
    stop('Please specify the correct type: "counts", "signature" or "density"')
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
  if (type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  } else {
    y.axis.labels <- round(y.axis.values, 0)

    # Write the mutation counts on top of graph
    for (i in 1 : 6) {
      j <- 23 + 32 * (i - 1)
      k <- 1 + 32 * (i - 1)
      text(bp[j], ymax * 7 / 8, labels = sum(counts[k : (32 * i)]),
           xpd = NA, cex = 0.8)
    }
  }
  text(-0.5, y.axis.values, labels = y.axis.labels,
       las = 1, adj = 1, xpd = NA, cex = cex)

  # Draw the x axis labels
  context.pos <- (bp[seq(1, 191, 2)] + bp[seq(2, 192, 2)]) / 2
  xlabel.1 <- c("A", "C", "G", "T")
  xlabel.2 <- rep(c("A", "C", "G", "T"), each = 4)
  text(context.pos, -ymax / 100, labels = rep(xlabel.1, 24), cex = 0.5,
       srt = 90, adj = 1, xpd = NA)
  text(context.pos, -ymax / 18, labels = rep(c("C", "T"), each = 48),
       cex = 0.5, srt = 90, adj = 1, xpd = NA)
  text(context.pos, -ymax / 10, labels = rep(xlabel.2, 6),
       cex = 0.5, srt = 90, adj = 1, xpd = NA)

  # Write the name of the sample
  text(1.5, ymax * 7 / 8, labels = id, adj = 0, cex = cex, font = 2)

  invisible(TRUE)
}

#' @rdname CatalogToPdf
#' @export
Cat192ToPdf <- function(catalog, name, id = colnames(catalog),
                        type = "counts", cex = 0.8, abundance = NULL) {
  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE)

  n <- ncol(catalog)
  graphics::par(mfrow = c(8, 1), mar = c(2, 4, 2, 2), oma = c(3, 2, 1, 1))

  # Do recycling of the function parameters if a vector
  # with length more than one is not specified by the user.
  if (n > 1 && length(type) == 1) {
    type <- rep(type, n)
  }

  for (i in 1 : n) {
    PlotCat192(catalog[, i, drop = FALSE],
               id = id[i], type = type[i],
               cex = cex, abundance = abundance)
  }
  invisible(grDevices::dev.off())
  invisible(TRUE)
}

#' @rdname PlotCatalog
#' @import graphics
#' @export
PlotCat192Strand <- function(catalog, id = colnames(catalog), type = "counts",
                             cex = 1, abundance = NULL) {
  stopifnot(dim(catalog) == c(192, 1))

  strand.col <- c('#394398',
                  '#e83020')

  # Sort data in plotting order
  counts <- catalog[.to.reorder.192.for.plotting, ]

  # Get the counts for each major mutation class
  counts.strand <- integer(12)
  for (i in 1 : 6){
    counts.strand[2 * i - 1] <-
      sum(counts[seq(32 * (i - 1) + 1, by = 2, length.out = 16)])
    counts.strand[2 * i] <-
      sum(counts[seq(32 * (i - 1) + 2, by = 2, length.out = 16)])
  }

  num.classes <- length(counts.strand)
  maj.class.names <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  cols <- rep(strand.col, num.classes / 2)

  if (type == "counts") {
    # Get ylim
    ymax <- max(counts.strand) * 1.3

    # Barplot: side by side
    mat <- matrix(counts.strand, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 5.5),
                  width = 0.3, xaxs = "i", yaxs = "i",
                  axes = FALSE, ann = FALSE, ylab = "counts",
                  border = NA, col = cols, xpd = NA)
  } else if (type == "signature") {
    # Calculate mutation signatures of each major mutation class
    sig <- counts.strand / sum(counts.strand)

    # Get ylim
    ymax <- ifelse(max(sig) * 1.3 > 1, 1, max(sig) * 1.3)

    # Barplot: side by side
    mat <- matrix(sig, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 5.5),
                  width = 0.3, xaxs = "i", yaxs = "i",
                  axes = FALSE, ann = FALSE, ylab = "proportion",
                  border = NA, col = cols, xpd = NA)
  } else if (type == "density") {
    stop('type = "density" not implemented')
  } else {
    stop('Please specify the correct type: "counts", "signature" or "density"')
  }

  # Draw y axis
  y.axis.values <- seq(0, ymax, length.out = 5)
  if (type != "counts") {
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

  # Draw the ID information on top of graph
  text(bp[5], ymax * 1.02, labels = id, xpd = NA,
       font = 2, cex = cex, adj = c(0, 0))

  invisible(TRUE)
}

#' @rdname CatalogToPdf
#' @export
Cat192StrandToPdf <- function(catalog, name, id = colnames(catalog),
                              type = "counts", cex = 1, abundance = NULL) {
  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE)

  n <- ncol(catalog)
  graphics::par(mfrow = c(4, 3), mar = c(2, 5, 2, 1), oma = c(2, 2, 2, 2))

  # Do recycling of the function parameters if a vector
  # with length more than one is not specified by the user.
  if (n > 1 && length(type) == 1) {
    type <- rep(type, n)
  }

  for (i in 1 : n) {
    PlotCat192Strand(catalog[, i, drop = FALSE],
                     id = id[i], type = type[i],
                     cex = cex, abundance = abundance)
  }
  invisible(grDevices::dev.off())

  invisible(TRUE)
}

#' @rdname PlotCatalog
#' @import graphics
#' @export
PlotCat1536 <- function(catalog, abundance, id = colnames(catalog)) {
  stopifnot(dim(catalog) == c(1536, 1))

  # Define the bases and their colors in plot
  bases <- c("A", "C", "G", "T")
  base.cols <- c("forestgreen", "dodgerblue2", "black", "red")

  # Define the theme color for plotting
  theme.col <- "darkgreen"
  colPal <- grDevices::colorRampPalette(c("white", theme.col))

  scale.col <- function(x, x.max) {
    idx <- round(x / x.max * 256)
    col <- colPal(256)[idx]
    return(col)
  }

  DrawImage <- function(counts, colors) {
    image(1 : 16, 1 : 16, matrix(counts, 16, 16),
          col = colors, axes = FALSE, ann = FALSE, bty = "n")

    # Draw the gridlines
    segments(rep(0.5, 5), c(0.5, 4.5, 8.5, 12.5, 16.5),
             rep(16.5, 5), c(0.5, 4.5, 8.5, 12.5, 16.5), xpd = T)
    segments(c(0.5, 4.5, 8.5, 12.5, 16.5), rep(0.5, 5),
             c(0.5, 4.5, 8.5, 12.5, 16.5), rep(16.5, 5), xpd = T)
  }

  DrawAxisY <- function() {
    text(0, 16 : 1, rep(bases, each = 4), col = rep(base.cols, each = 4),
         srt = 90, font = 2, xpd = T)
    text(-1, 16 : 1, rep(bases, 4), col = rep(base.cols, 4),
         srt = 90, font = 2, xpd = T)
    text(-2.5, 8.5, "Preceding bases", srt = 90, cex = 1.5, xpd = T)
    segments(rep(-1.3, 5), c(0.5, 4.5, 8.5, 12.5, 16.5), rep(0.5, 5),
             c(0.5, 4.5, 8.5, 12.5, 16.5), xpd = T)
  }

  DrawAxisX <- function() {
    text(1 : 16, 17, rep(bases, each = 4), col=rep(base.cols, each = 4),
         font = 2, xpd = T)
    text(1 : 16, 18, rep(bases, 4), col = rep(base.cols, 4),
         font = 2, xpd = T)
    segments(c(0.5, 4.5, 8.5, 12.5, 16.5), rep(16.5, 5),
             c(0.5, 4.5, 8.5, 12.5, 16.5), rep(18.3, 5), xpd = T)
  }

  mut.type <- rownames(catalog)

  # Calculate pentanucleotide sequence contexts, normalized by pentanucleotide
  # occurrence in the genome
  rates <- catalog
  for (i in 1 : 1536) {
    penta.names <- substr(mut.type[i], 1, 5)
    rates[i] <-
      catalog[i] * 1000000 / abundance[penta.names, ]
  }

  # Sort the rates matrix in plotting order
  rates <- as.data.frame(rates)
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
  old <- par(no.readonly = TRUE)
  par(mfrow = c(2, 3), oma = c(1, 1, 1, 1))

  for (i in 1 : n.types) {
    main.type <- main.types[i]
    sub.rates <- rates[main.mut.type == main.type, 1, drop = FALSE]

    # Do the color scaling
    max.col <- scale.col(max(sub.rates), max(rates))

    col.ref <- grDevices::colorRampPalette(c("white", max.col))(256)

    # Draw the 6 plots on page one by one
    if (i == 1) {
      par(mar = c(1, 5, 7.5, 1))
      DrawImage(sub.rates, col.ref)
      DrawAxisY()
      DrawAxisX()
      text(8.5, 19, main.type, cex = 1.5, xpd = T)
    }

    if (i == 2) {
      par(mar = c(1, 1, 7.5, 1))
      DrawImage(sub.rates, col.ref)
      DrawAxisX()
      text(8.5, 19, main.type, cex = 1.5, xpd = T)
      text(8.5, 20.5, id, cex = 1.5, xpd = T)
    }

    if (i == 3) {
      par(mar = c(1, 1, 7.5, 3))
      DrawImage(sub.rates, col.ref)
      DrawAxisX()
      text(8.5, 19, main.type, cex = 1.5, xpd = T)
      text(17.5, 17, '1bp 3\'', xpd = T, cex = 1)
      text(17.5, 18, '2bp 3\'', xpd = T, cex = 1)
    }

    if (i == 4) {
      par(mar = c(2.5, 5, 2, 1))
      DrawImage(sub.rates, col.ref)
      DrawAxisY()
      text(8.5, 17, main.type, cex = 1.5, xpd = T)
      text(-1, -0.7, '1bp 5\'', xpd = T, srt = 45, adj = 0, cex = 1)
      text(-2, -0.7, '2bp 5\'', xpd = T, srt = 45, adj = 0, cex = 1)
    }

    if (i == 5) {
      par(mar = c(2.5, 1, 2, 1))
      DrawImage(sub.rates, col.ref)
      text(8.5, 17, main.type, cex = 1.5, xpd = T)
    }

    if (i == 6) {
      par(mar = c(2.5, 1, 2, 3))
      DrawImage(sub.rates, col.ref)
      text(8.5, 17, main.type, cex = 1.5, xpd = T)
    }
  }
  on.exit(par(old), add = TRUE)

  invisible(TRUE)
}

#' @rdname CatalogToPdf
#' @export
Cat1536ToPdf <- function(catalog, name, id = colnames(catalog), abundance) {

  grDevices::cairo_pdf(name, width = 11.6929, height = 9.2677, onefile = TRUE)

  n <- ncol(catalog)

  for (i in 1 : n) {
    PlotCat1536(catalog[, i, drop = FALSE],
                id = id[i],
                abundance = abundance)
  }
  invisible(grDevices::dev.off())

  invisible(TRUE)
}

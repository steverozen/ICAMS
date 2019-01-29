#' @include plot_SNS_catalog.R
NULL

#' @rdname PlotCatalog
#' @import graphics
#' @export
PlotCatDNS78 <- function(catalog, id = colnames(catalog), type = "density",
                         abundance = NULL) {
  stopifnot(dim(catalog) == c(78, 1))
  stopifnot(rownames(catalog) == .catalog.row.order.DNS.78)

  dinuc.class.col <- RColorBrewer::brewer.pal(10, "Paired")
  cols <- rep(dinuc.class.col, c(9, 6, 9, 6, 9, 6, 6, 9, 9, 9))
  num.classes <- length(catalog)
  maj.class.names <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")

  if (type == "density") {
    # Calculate rate of mutations per million nucleotides for the catalog
    rate <- double(78)
    for (i in 1 : 78) {
      rate[i] <-
        catalog[i] * 1000000 / abundance[substr(rownames(catalog)[i], 1, 2), ]
    }

    # Get ylim
    ymax <- ifelse(max(rate) * 1.3 > 1, 1, max(rate) * 1.3)

    # Barplot
    bp <- barplot(rate, ylim = c(0, ymax), xaxt = "n", yaxt = "n", xaxs = "i",
                  lwd = 3, space = 1.35, border = NA, col = cols,
                  xpd = NA, ylab = "mut/million")

    # Write the mutation counts on top of graph
    for (i in 1 : 10) {
      j <- c(7, 13, 20, 28, 35, 43, 49, 56, 65, 74)
      name <- substr(rownames(catalog), 1, 2)
      text(bp[j[i]], ymax * 0.9, xpd = NA, cex = 0.8,
           labels = sum(catalog[name == maj.class.names[i], ]))
    }
  } else if (type == "counts") {
    # Get ylim
    ymax <- max(catalog[, 1]) * 1.3

    # Barplot
    bp <- barplot(catalog[, 1], xaxt = "n", yaxt = "n", ylim = c(0, ymax),
                  xaxs = "i", lwd = 3, space = 1.35, border = NA, xaxs = "i",
                  col = cols, ylab = "counts")

    # Write the mutation counts on top of graph
    for (i in 1 : 10) {
      j <- c(7, 13, 20, 28, 35, 43, 49, 56, 65, 74)
      name <- substr(rownames(catalog), 1, 2)
      text(bp[j[i]], ymax * 0.9, xpd = NA, cex = 0.8,
           labels = sum(catalog[name == maj.class.names[i], ]))
    }
  } else if (type == "signature") {
    # Calculate mutation signatures of the input catalog
    sig <- catalog / sum(catalog)

    # Get ylim
    ymax <- ifelse(max(sig) * 1.3 > 1, 1, max(sig) * 1.3)

    # Barplot
    bp <- barplot(sig[, 1], xaxt = "n", yaxt = "n", ylim = c(0, ymax),
                  lwd = 3, space = 1.35, border = NA, xaxs = "i",
                  col = cols, ylab = "proportion")
  } else {
    stop('Please specify the correct type: "density", "counts" or "signature"')
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

  # Draw the ID information on top of graph
  text(1.5, ymax * 7 / 8, labels = id, adj = 0, cex = 0.8, font = 2)

  # Draw y axis
  y.axis.values <- seq(0, ymax, ymax / 4)
  if (type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  } else {
    y.axis.labels <- round(y.axis.values, 0)
  }
  text(0.35, y.axis.values, labels = y.axis.labels,
       las = 1, adj = 1, xpd = NA, cex = 0.75)

  # Draw x axis labels
  text(bp, -ymax / 80, labels = substr(rownames(catalog), 4, 4),
       cex = 0.5, srt = 90, adj = 1, xpd = NA)
  text(bp, -ymax / 15, labels = substr(rownames(catalog), 3, 3),
       cex = 0.5, srt = 90, adj = 1, xpd = NA)

  invisible(TRUE)
}

#' @rdname CatalogToPdf
#' @export
CatDNS78ToPdf <-
  function(catalog, name, id = colnames(catalog),
           type = "density", abundance = NULL) {
    # Setting the width and length for A4 size plotting
    grDevices::cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE)

    n <- ncol(catalog)
    graphics::par(mfrow = c(8, 1), mar = c(2, 4, 2, 2), oma = c(3, 3, 2, 2))

    # Do recycling of the function parameters if a vector
    # with length more than one is not specified by the user.
    if (n > 1 && length(type) == 1) {
      type <- rep(type, n)
    }

    for (i in 1 : n) {
      PlotCatDNS78(catalog[, i, drop = FALSE],
                   id = id[i],
                   type = type[i],
                   abundance = abundance)
    }
    invisible(grDevices::dev.off())

    invisible(TRUE)
  }

#' @rdname PlotCatalog
#' @import graphics
#' @export
PlotCatDNS144 <- function(catalog, id = colnames(catalog), type = "counts",
                          cex = 1, abundance = NULL) {
  stopifnot(dim(catalog) == c(144, 1))

  strand.col <- c('#394398',
                  '#e83020')

  # Sort data in plotting order
  counts <- catalog[.to.reorder.144.for.plotting, ]

  # Get the counts for each major mutation class
  counts.strand <- integer(20)
  for (i in 1 : 10){
    idx <- c(0, 18, 24, 42, 48, 66, 72, 78, 96, 114, 132)
    counts.strand[2 * i - 1] <-
      sum(counts[seq(idx[i] + 1, idx[i + 1], by = 2)])
    counts.strand[2 * i] <-
      sum(counts[seq(idx[i] + 2, idx[i + 1], by = 2)])
  }

  num.classes <- length(counts.strand)
  maj.class.names <-
    c("AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT")
  cols <- rep(strand.col, num.classes / 2)

  if (type == "counts") {
    # Get ylim
    ymax <- max(counts.strand) * 1.3

    # Barplot: side by side
    mat <- matrix(counts.strand, nrow = 2, ncol = num.classes / 2)
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 9),
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
    bp <- barplot(mat, beside = TRUE, ylim = c(0, ymax), xlim = c(0, 9),
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

  # Draw the ID information on top of graph
  text(bp[8], ymax, labels = id, xpd = NA,
       font = 2, cex = cex, adj = c(0, 0))

  invisible(TRUE)
}

#' @rdname CatalogToPdf
#' @export
CatDNS144ToPdf <- function(catalog, name, id = colnames(catalog),
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
    PlotCatDNS144(catalog[, i, drop = FALSE],
                  id = id[i], type = type[i],
                  cex = cex, abundance = abundance)
  }
  invisible(grDevices::dev.off())

  invisible(TRUE)
}

#' @rdname PlotCatalog
#' @import graphics
#' @export
PlotCatQUAD136 <- function(catalog, id = colnames(catalog),
                           type = "density", abundance = NULL) {
  stopifnot(dim(catalog) == c(136, 1))

  # Specify the lay out of the plotting
  invisible(layout(matrix(c(7, 8, 9, 10, 4, 5, 6, 11, 1, 2 , 3, 11), 3, 4,
                          byrow = TRUE)))

  # Define the bases and their colors in plot
  base <- c("A", "C", "G", "T")
  base.cols <- c("forestgreen", "dodgerblue2", "black", "red")

  ref.order <- c("AC", "AT", "GC", "CC", "CG", "CT", "TA", "TC", "TG", "TT")
  mut.type <- paste(ref.order, "NN", sep = ">")

  # Calculate the occurrences of each mutation type for plotting
  counts <- matrix(0, nrow = 160, ncol = 1)
  rownames(counts) <- .order.for.QUAD136.plotting
  for (i in 1:160){
    if (.order.for.QUAD136.plotting[i] %in% rownames(catalog)) {
      counts[i] <- catalog[.order.for.QUAD136.plotting[i], ]
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

  if (type == "density") {
    # Calculate tetranucleotide sequence contexts, normalized by tetranucleotide
    # occurrence in the genome
    rates <- matrix(0, nrow = 160, ncol = 1)
    rownames(rates) <- .order.for.QUAD136.plotting
    for (i in 1:160){
      if (.order.for.QUAD136.plotting[i] %in% rownames(catalog)) {
        rates[i] <-
          catalog[.order.for.QUAD136.plotting[i], ] /
          abundance[.order.for.QUAD136.plotting[i], ]
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
  }

  DrawImage <- function(mat) {
    value <- as.numeric(unlist(mat))
    maximum <- max(value[!is.na(value)])
    if (maximum == 0) {
      col.ref <- "white"
    } else {
      col.ref <- "palegreen3"
    }

    image(1:4, 1:4, mat,
          col = grDevices::colorRampPalette(c("white", col.ref))(16),
          asp = 1, axes = FALSE, ann = FALSE)

    # Make the background of the plot grey
    rect(0.5, 0.5, 4.5, 4.5 , col = "grey")

    # Plot the image again
    image(1:4, 1:4, mat,
          col = grDevices::colorRampPalette(c("white", col.ref))(16),
          asp = 1, axes = FALSE, ann = FALSE, add = TRUE)
  }

  for (i in 1:10){
    par(mar = c(1, 1, 4.5, 1))

    if (type == "density") {
      DrawImage(matrix(rates[(16 * (i - 1) + 1) : (16 * i)], 4, 4))
    } else if (type == "counts") {
      DrawImage(matrix(counts[(16 * (i - 1) + 1) : (16 * i)], 4, 4))
    } else {
      stop('Please specify the correct type: "density" or "counts"')
    }

    # Draw the mutation type and number of occurrences on top of image
    text(2, 5.3, mut.type[i], font = 2, xpd = NA)
    text(3.2, 5.3, paste0("(", counts.per.class[ref.order[i], ], ")"), font = 2, xpd = NA)

    # Draw a box surrounding the image
    segments(c(0.5, 0.5), c(0.5, 4.5), c(4.5, 4.5), c(0.5, 4.5), xpd = NA)
    segments(c(0.5, 4.5), c(0.5, 0.5), c(0.5, 4.5), c(4.5, 4.5), xpd = NA)

    # Draw the base information of the plot
    text(rep(0.2, 4), c(4, 3, 2, 1), base, col = base.cols, xpd = NA)
    text(seq(4), rep(4.8, 4), base, col = base.cols, xpd = NA)

    # Draw the ID information of the sample
    if (i == 8) {
      mtext(id, at = 5, line =3)
    }
  }

  # Add in additional information
  plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
  text(x = 0.5, y = 0.9, "Maxima per class", cex = 1.6)
  ref <- c("TA", "TC", "TG", "TT", "CC", "CG", "CT", "AC", "AT", "GC")

  if (type == "density") {
    text(x = 0.5, y = 0.8, "(mut/million)", cex = 1.2)
    maxima <- numeric(0)
    for (i in 1:10) {
      maxima[i] <- max.rate.per.class[ref[i], ]
      names(maxima)[i] <- ref[i]
    }
  } else if (type == "counts") {
    text(x = 0.5, y = 0.8, "(counts)", cex = 1.2)
    maxima <- numeric(0)
    for (i in 1:10) {
      maxima[i] <- max.count.per.class[ref[i], ]
      names(maxima)[i] <- ref[i]
    }
  } else {
    stop('Please specify the correct type: "density" or "counts"')
  }
  text(rep(0, 5), seq(0.7, 0.3, length.out = 5),
       paste(ref[1:5], maxima[1:5], sep = " = "), adj = 0, cex = 1.2)
  text(rep(0.6, 5), seq(0.7, 0.3, length.out = 5),
       paste(ref[6:10], maxima[6:10], sep = " = "), adj = 0, cex = 1.2)
}

#' @rdname CatalogToPdf
#' @export
CatQUAD136ToPdf <- function(catalog, name, id = colnames(catalog),
                            type = "density", abundance = NULL) {
  stopifnot(nrow(catalog) == 136)
  n <- ncol(catalog)

  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE)
  par(oma = c(2, 2, 2, 2))

  # Do recycling of the function parameters if a vector
  # with length more than one is not specified by the user.
  if (n > 1 && length(type) == 1) {
    type <- rep(type, n)
  }

  # Specify the lay out of the plotting
  invisible(layout(matrix(c(7, 8, 9, 10, 4, 5, 6, 11, 1, 2, 3, 11, 18,
                            19, 20, 21, 15, 16, 17, 22, 12, 13, 14, 22),
                          6, 4,byrow = TRUE)))

  # Define the bases and their colors in plot
  base <- c("A", "C", "G", "T")
  base.cols <- c("forestgreen", "dodgerblue2", "black", "red")

  ref.order <- c("AC", "AT", "GC", "CC", "CG", "CT", "TA", "TC", "TG", "TT")
  mut.type <- paste(ref.order, "NN", sep = ">")

  for (i in 1:n) {
    cat <- catalog[, i, drop = FALSE]

    # Calculate the occurrences of each mutation type for plotting
    counts <- matrix(0, nrow = 160, ncol = 1)
    rownames(counts) <- .order.for.QUAD136.plotting
    for (j in 1:160){
      if (.order.for.QUAD136.plotting[j] %in% rownames(cat)) {
        counts[j] <- cat[.order.for.QUAD136.plotting[j], ]
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

    if (type[i] == "density") {
      # Calculate tetranucleotide sequence contexts, normalized by tetranucleotide
      # occurrence in the genome
      rates <- matrix(0, nrow = 160, ncol = 1)
      rownames(rates) <- .order.for.QUAD136.plotting
      for (j in 1:160){
        if (.order.for.QUAD136.plotting[j] %in% rownames(cat)) {
          rates[j] <-
            cat[.order.for.QUAD136.plotting[j], ] /
            abundance[.order.for.QUAD136.plotting[j], ]
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
        col.ref <- "palegreen3"
      }

      image(1:4, 1:4, mat,
            col = grDevices::colorRampPalette(c("white", col.ref))(16),
            asp = 1, axes = FALSE, ann = FALSE)

      # Make the background of the plot grey
      rect(0.5, 0.5, 4.5, 4.5 , col = "grey")

      # Plot the image again
      image(1:4, 1:4, mat,
            col = grDevices::colorRampPalette(c("white", col.ref))(16),
            asp = 1, axes = FALSE, ann = FALSE, add = TRUE)
    }

    for (j in 1:10) {
      par(mar = c(1, 1, 4.5, 1))
      if (type[i] == "density") {
        DrawImage(matrix(rates[(16 * (j - 1) + 1) : (16 * j)], 4, 4))
      } else if (type[i] == "counts") {
        DrawImage(matrix(counts[(16 * (j - 1) + 1) : (16 * j)], 4, 4))
      } else {
        stop('Please specify the correct type: "density" or "counts"')
      }

      # Draw the mutation type and number of occurrences on top of image
      text(2, 5.3, mut.type[j], font = 2, xpd = NA)
      text(3.2, 5.3, paste0("(", counts.per.class[ref.order[j], ], ")"), font = 2, xpd = NA)

      # Draw a box surrounding the image
      segments(c(0.5, 0.5), c(0.5, 4.5), c(4.5, 4.5), c(0.5, 4.5), xpd = NA)
      segments(c(0.5, 4.5), c(0.5, 0.5), c(0.5, 4.5), c(4.5, 4.5), xpd = NA)

      # Draw the base information of the plot
      text(rep(0.2, 4), c(4, 3, 2, 1), base, col = base.cols, xpd = NA)
      text(seq(4), rep(4.8, 4), base, col = base.cols, xpd = NA)

      # Draw the ID information of the sample
      if (j == 8) {
        mtext(id[i], at = 5, line = 3)
      }
    }

    # Add in additional information
    plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
    text(x = 0.5, y = 0.9, "Maxima per class", cex = 1.6)
    ref <- c("TA", "TC", "TG", "TT", "CC", "CG", "CT", "AC", "AT", "GC")

    if (type[i] == "density") {
      text(x = 0.5, y = 0.8, "(mut/million)", cex = 1.2)
      maxima <- numeric(0)
      for (j in 1:10) {
        maxima[j] <- max.rate.per.class[ref[j], ]
        names(maxima)[j] <- ref[j]
      }
    } else if (type[i] == "counts") {
      text(x = 0.5, y = 0.8, "(counts)", cex = 1.2)
      maxima <- numeric(0)
      for (j in 1:10) {
        maxima[j] <- max.count.per.class[ref[j], ]
        names(maxima)[j] <- ref[j]
      }
    } else {
      stop('Please specify the correct type: "density" or "counts"')
    }
    text(rep(0, 5), seq(0.7, 0.3, length.out = 5),
         paste(ref[1:5], maxima[1:5], sep = " = "), adj = 0, cex = 1.2)
    text(rep(0.6, 5), seq(0.7, 0.3, length.out = 5),
         paste(ref[6:10], maxima[6:10], sep = " = "), adj = 0, cex = 1.2)
  }
  invisible(grDevices::dev.off())
  invisible(TRUE)
}

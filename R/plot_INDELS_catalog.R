#' @include plot_SNS_catalog.R
NULL

#' @rdname PlotCatalog
#' @import graphics
#' @export
PlotCatID <- function(catalog, id = colnames(catalog), type = "counts"){
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

  if (type == "counts") {
    # Get ylim
    ymax <- max(catalog) * 1.3

    # Barplot
    bp <- barplot(catalog[, 1], ylim = c(0, ymax), axes = FALSE, xaxt = "n",
                  lwd = 3, space = 1.35, border = NA, col = cols, xpd = NA,
                  xaxs = "i", ylab = "counts")
  } else if (type == "signature") {
    # Calculate mutation signatures of the input catalog
    sig <- catalog / sum(catalog)

    # Get ylim
    ymax <- ifelse(max(sig) * 1.3 > 1, 1, max(sig) * 1.3)

    # Barplot
    bp <- barplot(sig[, 1], ylim = c(0, ymax), axes = FALSE, xaxt = "n",
                  lwd = 3, space = 1.35, border = NA, col = cols, xpd = NA,
                  xaxs = "i", ylab = "proportion")
  } else {
    stop('Please specify the correct type: "counts" or "signature"')
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
                 (x.left[c(5, 9)] + x.right[c(8, 12)] - 12) / 2,
                 (x.left[13] + x.right[length(x.left)]) / 2)
  category.lab <- c(rep(c("C", "T"), 2), rep(c("2", "3", "4", "5+"), 3))
  category.col <- c(rep(c("black", "white"), 2),
                    rep(c("black", "black", "black", "white"), 3))

  # Draw lines above each class
  rect(xleft = x.left, ymax * 1.01, xright = x.right, ymax * 1.09,
       col = indel.class.col, border = NA, xpd = NA)
  text((x.left + x.right) / 2, ymax * 1.05, labels = category.lab,
       cex = 0.55, col = category.col, xpd = NA)

  # Draw mutation class labels at the top of the figure
  text(class.pos, ymax * 1.23, labels = maj.class.names, cex = 0.75, xpd = NA)

  # Draw the ID information of the sample
  text(1.5, ymax * 7 / 8, labels = id, adj = 0, cex = 0.85, font = 2)

  # Draw y axis
  y.axis.values <- seq(0, ymax, ymax / 4)
  if (type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  } else {
    y.axis.labels <- round(y.axis.values, 0)
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

  invisible(TRUE)
}

#' @rdname CatalogToPdf
#' @export
CatIDToPdf <-
  function(catalog, name, id = colnames(catalog), type = "counts") {
    # Setting the width and length for A4 size plotting
    grDevices::cairo_pdf(name, width = 8.2677, height = 11.6929, onefile = TRUE)

    n <- ncol(catalog)
    graphics::par(mfrow = c(8, 1), mar = c(3, 4, 2.5, 2), oma = c(3, 3, 2, 2))

    # Do recycling of the function parameters if a vector
    # with length more than one is not specified by the user.
    if (n > 1 && length(type) == 1) {
      type <- rep(type, n)
    }

    for (i in 1 : n) {
      PlotCatID(catalog[, i, drop = FALSE],
                id = id[i],
                type = type[i])
    }
    invisible(grDevices::dev.off())

    invisible(TRUE)
  }

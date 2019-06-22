#' Create position probability matrix (PPM) for *one* sample from
#' a Variant Call Format (VCF) file.
#'
#' @param vcf One in-memory data frame of pure SBS mutations -- no DBS or 3+BS
#'   mutations.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param seq.context.width The number of preceding and following bases to be
#'   extracted around the mutated position from \code{ref.genome}.
#'
#' @importFrom utils tail
#'
#' @return A position probability matrix (PPM). 
#'
#' @keywords internal
CreateOnePPMFromSBSVCF <- function(vcf, ref.genome, seq.context.width) {
  stopifnot(nrow(vcf) != 0)
  vcf <- AddSeqContext(vcf, ref.genome = ref.genome,
                       seq.context.width = seq.context.width)
  dt <- data.table(vcf)

  # Map the sequence context column in dt to strand-agnostic category
  idx <- substr(dt[[tail(names(dt), 1)]], seq.context.width + 1,
                seq.context.width + 1) %in% c("A", "G")
  dt[[tail(names(dt), 1)]][idx] <- revc(dt[[tail(names(dt), 1)]][idx])

  # Create the position probability matrix (PPM)
  GetPPM <- function(idx, seq.context) {
    base <- substr(seq.context, idx, idx)
    count <- table(factor(base, levels = c("A", "C", "G", "T")))
    return(as.matrix(count / sum(count)))
  }
  mat <- sapply(1:nchar(dt[[tail(names(dt), 1)]][1]), FUN = GetPPM,
                seq.context = dt[[tail(names(dt), 1)]])
  rownames(mat) <- c("A", "C", "G", "T")
  colnames(mat) <- c(paste0(seq.context.width:1, "bp 5'"), "target",
                     paste0(1:seq.context.width, "bp 3'"))
  return(mat)
}

#' Create position probability matrices (PPM) from a list of SBS vcfs
#'
#' @param list.of.SBS.vcfs List of in-memory data frames of pure SBS mutations
#'   -- no DBS or 3+BS mutations.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param seq.context.width The number of preceding and following bases to be
#'   extracted around the mutated position from \code{ref.genome}.
#'
#' @return A list of position probability matrices (PPM).
#'
#' @keywords internal
CreatePPMFromSBSVCFs <-
  function(list.of.SBS.vcfs, ref.genome, seq.context.width) {
    list.of.PPM <- lapply(list.of.SBS.vcfs, FUN = CreateOnePPMFromSBSVCF,
                          ref.genome = ref.genome,
                          seq.context.width = seq.context.width)
    return(list.of.PPM)
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
  graphics::par(mfrow = c(4, 2), mar = c(4, 5.5, 2, 1), oma = c(1, 1, 2, 1))

  for (i in 1:n) {
    PlotPPM(list.of.ppm[[i]], title = titles[i])
  }
  invisible(grDevices::dev.off())
  invisible(TRUE)
}

#' Create position weight matrix (PWM) for *one* sample from
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
#' @return A position weight matrix (PWM) with each row representing the
#'   position weight of four bases.
#'
#' @keywords internal
CreateOnePWMFromSBSVCF <- function(vcf, ref.genome, seq.context.width) {
  stopifnot(nrow(vcf) != 0)
  vcf <- AddSeqContext(vcf, ref.genome = ref.genome,
                       seq.context.width = seq.context.width)
  dt <- data.table(vcf)

  # Map the sequence context column in dt to strand-agnostic category
  idx <- substr(dt[[tail(names(dt), 1)]], seq.context.width + 1,
                seq.context.width + 1) %in% c("A", "G")
  dt[[tail(names(dt), 1)]][idx] <- revc(dt[[tail(names(dt), 1)]][idx])

  # Create the position weight matrix (PWM)
  GetPWM <- function(idx, seq.context) {
    base <- substr(seq.context, idx, idx)
    count <- table(factor(base, levels = c("A", "C", "G", "T")))
    return(as.matrix(count / sum(count)))
  }
  mat <- sapply(1:nchar(dt[[tail(names(dt), 1)]][1]), FUN = GetPWM,
                seq.context = dt[[tail(names(dt), 1)]])
  rownames(mat) <- c("A", "C", "G", "T")
  colnames(mat) <- c(paste0(seq.context.width:1, "bp 5'"), "target",
                     paste0(1:seq.context.width, "bp 3'"))
  return(t(mat))
}

#' Create position weight matrices (PWM) from a list of SBS vcfs
#'
#' @param list.of.SBS.vcfs List of in-memory data frames of pure SBS mutations
#'   -- no DBS or 3+BS mutations. The list names will be the main title on top
#'   of each position weight matrices (PWM) plot.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param ref.context.width Width of the reference context to be used for
#'   generating the position weight matrix (PWM).
#'
#' @return A list of position weight matrices (PWM).
#'
#' @keywords internal
CreatePWMFromSBSVCFs <-
  function(list.of.SBS.vcfs, ref.genome, ref.context.width) {
    list.of.SBS.vcfs <- lapply(list.of.SBS.vcfs, FUN = AddSeqContext,
                               ref.genome = ref.genome)
    list.of.PWM <- lapply(list.of.SBS.vcfs, FUN = CreateOnePWMFromSBSVCF,
                          ref.context.width = ref.context.width)
    return(list.of.PWM)
  }

#' Plot position weight matrix (PWM) for *one* sample from a Variant Call Format
#' (VCF) file.
#'
#' @param pwm A position weight matrix (PWM) for *one* sample with each row
#'   representing the position weight of four bases.
#'
#' @param title The main title of the plot.
#'
#' @return \code{invisible(TRUE)}
#'
#' @keywords internal
PlotPWM <- function(pwm, title) {
  par(mar = c(5, 5, 1, 1))
  x <- seq(1:nrow(pwm))
  plot(x, y = pwm[, "A"], xaxt = "n", xlab = "",
       ylab = "frequency", ylim = c(0, 1), pch = 20, type = "b",
       lwd = 2, col = "darkgreen")
  lines(x, pwm[, "C"], col = "blue", type = "b", pch = 20, lwd = 2)
  lines(x, pwm[, "G"], col = "black", type = "b", pch = 20, lwd = 2)
  lines(x, pwm[, "T"], col = "red", type = "b", pch = 20, lwd = 2)
  abline(h = 0.25, col = "grey50")
  text(x, y = -0.1, labels = rownames(pwm), adj = 1, srt = 90, xpd = NA)
  legend("topright", pch = 16, ncol = 4, bty = "n",
         legend = c("A", "C", "G", "T"),
         col = c("darkgreen", "blue", "black", "red"))
  mtext(text = title, side = 3, line = 0.5)
  invisible(TRUE)
}

#' Plot position weight matrices (PWM) to a PDF file
#'
#' @param list.of.pwm A list of position weight matrices (PWM)
#'
#' @param file The name of the PDF file to be produced.
#'
#' @param titles A vector of titles on top of each PWM plot.
#'
#' @return \code{invisible(TRUE)}
#'
#' @keywords internal
PlotPWMToPdf <- function(list.of.pwm, file, titles = names(list.of.pwm)) {
  # Setting the width and length for A4 size plotting
  grDevices::cairo_pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)

  n <- length(list.of.pwm)
  graphics::par(mfrow = c(2, 1), mar = c(4, 5.5, 2, 1), oma = c(1, 1, 2, 1))

  for (i in 1:n) {
    PlotPWM(list.of.pwm[[i]], title = titles[i])
  }
  invisible(grDevices::dev.off())
  invisible(TRUE)
}

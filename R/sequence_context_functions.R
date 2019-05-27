#' Create position weight matrix (PWM) for *one* sample from
#' a Variant Call Format (VCF) file.
#'
#' @param vcf An in-memory VCF file annotated by the AddSequence function. It
#'   must *not* contain indels and must *not* contain DBS (double base
#'   substitutions), or triplet base substitutions etc., even if encoded as
#'   neighboring SBS.
#'
#' @param ref.context.width Width of the reference context to be used for
#'   generating the position weight matrix (PWM).
#'
#' @return A position weight matrix (PWM) with each row representing the
#'   position weight of four bases.
#'
#' @keywords internal
CreateOnePWMFromSBSVCF <- function(vcf, ref.context.width) {
  stopifnot(nrow(vcf) != 0)
  stopifnot("seq.21context" %in% colnames(vcf))
  dt <- data.table(vcf)

  # Map the seq.21context column in dt to strand-agnostic category
  idx <- substr(dt$seq.21context, 11, 11) %in% c("A", "G")
  dt$seq.21context[idx] <- revc(dt$seq.21context[idx])

  # Get the desired reference context based on ref.context.width
  k <- (ref.context.width - 1) / 2
  dt[, ref.context := substr(seq.21context, 11 - k, 11 + k)]

  # Create the position weight matrix (PWM)
  GetPWM <- function(idx, ref.context) {
    base <- substr(ref.context, idx, idx)
    count <- table(factor(base, levels = c("A", "C", "G", "T")))
    return(as.matrix(count / sum(count)))
  }
  mat <- sapply(1:ref.context.width, FUN = GetPWM,
                ref.context = dt$ref.context)
  rownames(mat) <- c("A", "C", "G", "T")
  colnames(mat) <- c(paste0(k:1, "bp 5'"), "target", paste0(1:k, "bp 3'"))
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
    list.of.SBS.vcfs <- lapply(list.of.SBS.vcfs, FUN = AddSequence,
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

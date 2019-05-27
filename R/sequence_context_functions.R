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


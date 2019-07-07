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

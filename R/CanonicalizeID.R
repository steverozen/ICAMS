#' @title Determine the mutation types of insertions and deletions.
#'
#' @param context A vector of ample surrounding
#'   sequence on each side the variants
#'
#' @param ref Vector of reference alleles; this includes
#' one un-altered base at the start of the reference allele
#' e.g. for a deletion of a single T this might be "AT", in
#' which case the \code{alt} allele woult be "A"
#'
#' @param alt Vector of alternative alleles.
#'
#' @param pos Vector of the positions of the insertions and deletions in
#'  \code{context}. This is the position of the unaltered allele shared
#' by \code{ref} and \code{alt}. So in 1-based indexing, it is the
#' position just before the insertion or the deletion.
#'
#' @return A vector of strings that are the canonical representations
#'  of the given insertions and deletions.
#'
#' @importFrom utils head
#'
#' @keywords internal
CanonicalizeID <- function(context, ref, alt, pos) {
  ret <- mapply(
    justify_and_categorize_1_indel,
    context,
    ref,
    alt,
    orig_pos = pos,
    explain_indels = FALSE
  )
  return(ret)
}

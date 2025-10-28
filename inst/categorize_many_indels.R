# categorize_indel_mutations

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
#' @return A data frame parallel to the input vectors. The
#' data frame has the columns COSMIC_83, ins_or_del, previous_char, repeat_seq, repeat_count, post_char, mh_seq
#'
#' @importFrom utils head
#'
#' @keywords internal
categorize_many_indels <- function(vcf) {
  context = vcf$seq.context
  ref = vcf$REF
  alt = vcf$ALT
  pos = vcf$seq.context.width + 1
  if (all(substr(ref, 1, 1) == substr(alt, 1, 1))) {
    ref <- substr(ref, 2, nchar(ref))
    alt <- substr(alt, 2, nchar(alt))
    pos = pos + 1
    # + 1 because ref and alt shared the first, shared character,
    # e.g. ref = "AG", alt = "A" deletion of G
    # e.g. ref = "T", alt = "TCC" -- insertion of CC
  } else {
    stopifnot(ref != "" | alt != "")
  }

  ret <- mapply(xCanonicalize1ID, context, ref, alt, pos, 0, regress = FALSE)

  return(ret)
}

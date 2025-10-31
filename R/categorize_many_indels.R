#' @title Determine the mutation types of insertions and deletions.
#'
#' @param vcf A dataframe represenging a VCF, which must contain the
#' following columns
#'
#' seq.context: ample surrounding
#'   sequence on each side the variants
#'
#' ref The reference alleles as strings; this includes
#' one un-altered base at the start of the reference allele
#' e.g. for a deletion of a single T this might be "AT", in
#' which case the \code{alt} allele woult be "A"
#'
#' alt The alternative alleles as strings
#'
#' pos The positions of the insertions and deletions in
#'  \code{context}. This is the position of the unaltered allele shared
#' by \code{ref} and \code{alt}. So in 1-based indexing, it is the
#' position just before the insertion or the deletion.
#'
#' @param explain_indels Generate verbose stdout messages regarding
#' categorization process
#'
#' @return A data frame parallel to the input vectors. The
#' data frame has the columns COSMIC_83, ins_or_del, previous_char, repeat_seq, repeat_count, post_char, mh_seq
#'
#' @md
#'
#' @keywords internal
categorize_many_indels <- function(vcf, explain_indels = FALSE) {
  context = vcf[, "seq.context"] # Make sure there's an error if there is no seq.context column
  ref = vcf$REF
  alt = vcf$ALT
  pos = vcf$seq.context.width + 1

  ret <- mapply(
    justify_and_categorize_1_indel,
    context = context,
    orig_ref = ref,
    orig_alt = alt,
    orig_pos = pos,
    explain_indels = explain_indels,
    SIMPLIFY = FALSE
  )

  cnames = colnames(vcf)
  if ("CHROM" %in% cnames && "POS" %in% cnames) {
    names(ret) = paste(vcf$CHROM, vcf$POS, sep = "_")
  } else {
    names(ret) = ret$COSMIC_83
  }

  return(ret)
}

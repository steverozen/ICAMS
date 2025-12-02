#' @title Given a single insertion or deletion in context categorize it.
#'
#' @param context Ample surrounding sequence on each side of the insertion or
#'   deletion.
#'
#' @param orig_ref The reference allele as a single string.
#'
#' @param orig_alt The alternative allele as a single string
#'
#' @param orig_pos The position of the inserted or deleted sequence in
#' \code{context}. If this is an insertion, it is the positon after the position
#' of the insertion of \code{alt}. If this is a deletion, it is the position in
#' context at which \code{ref} starts.
#'
#' @return A string that is the canonical representation of the type of the
#'  given insertion or deletion. Return \code{NA} and raise a warning if there
#'  is an un-normalized representation of the deletion of a repeat unit. See
#'  \code{FindDelMH} for details. (This seems to be very rare.)
#'
#' @keywords internal
#'
x_categorize_1_indel = function(
  context,
  orig_ref,
  orig_alt,
  orig_pos
) {
  requireNamespace("stringi")

  regress = (Sys.getenv("ICAMS_DO_NOT_REGRESS") == '')

  if (orig_ref == "" || orig_alt == "") {
    stop(
      "cannot remove common prefix if orig_ref or orig_alt is the empty string"
    )
  }
  stopifnot(substr(orig_ref, 1, 1) == substr(orig_alt, 1, 1))
  ref <- substr(orig_ref, 2, nchar(orig_ref))
  alt <- substr(orig_alt, 2, nchar(orig_alt))
  pos = orig_pos + 1

  if (nchar(alt) < nchar(ref)) {
    retval = categorize_1_justified_indel(
      context = context,
      ins_or_del = "d",
      ins_or_del_seq = ref,
      pos = pos
    )

    if (regress) {
      prev_ret = Canonicalize1Del(context, ref, pos) # pos is the start of the deletion

      if (is.na(prev_ret) || prev_ret != retval$COSMIC_83) {
        message("\n\nDELETION difference:")
        message("old = ", prev_ret)
        message("new = ", retval$COSMIC_83, " ref = ", ref)
      }
    }
  } else if (nchar(alt) > nchar(ref)) {
    retval = categorize_1_justified_indel(
      context = context, # Not sure if this is correct, of if need to do the insertion
      ins_or_del = "i",
      ins_or_del_seq = alt,
      pos = pos
    )

    if (regress) {
      prev_ret = Canonicalize1INS(context, alt, pos - 1) # the insertion occurs immediately after pos - 1

      if (prev_ret != retval$COSMIC_83) {
        message("\nINSERTION difference:")
        message("old = ", prev_ret)
        message("new = ", retval$COSMIC_83, "\nalt = ", alt)
      }
    }
  } else {
    message("Non-insertion / non-deletion found: ", ref, " ", alt, " ", context)
    return(indel_all_na_return())
  }
  if (regress) {
    retval$prev_COSMIC_83 = prev_ret
  }
  return(retval)
}

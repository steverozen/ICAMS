#' @title Given a single insertion or deletion in context categorize it.
#'
#' @param context Ample surrounding
#'   sequence on each side of the insertion or deletion.
#'
#' @param ref The reference allele as a single string.
#'
#' @param alt The alternative allele as a single string
#'
#' @param pos The position of the inserted or deleted sequence in \code{context}.
#' If this is an insertion, it is the positon after the position of the insertion of \code{alt}.
#' If this is a deletion, it is the position in context at which \code{ref} starts.
#'
#' @param trace If > 0, then generate messages tracing
#' how the computation is carried out.
#
#' @return A string that is the canonical representation
#'  of the type of the given
#'  insertion or deletion.
#'  Return \code{NA}
#'  and raise a warning if
#'  there is an un-normalized representation of
#'  the deletion of a repeat unit.
#'  See \code{FindDelMH} for details.
#'  (This seems to be very rare.)
#'
#' @keywords internal
xCanonicalize1ID <- function(
  context,
  ref,
  alt,
  pos,
  trace = 0,
  regress = TRUE
) {
  if (trace > 0) {
    message("Canonicalize1ID(", context, ",", ref, ",", alt, ",", pos, "\n")
  }
  if (nchar(alt) < nchar(ref)) {
    # We have a deletion.
    # pos is the 1-based start of the deletion in context.

    tmp_long = context
    end_pos = pos + nchar(ref) - 1
    tmp_short = stringi::stri_sub_replace(
      context,
      from = pos,
      to = end_pos,
      replacement = ""
    )
    tmp = justify_indel(tmp_long, tmp_short, pos, ref)
    new_ret = ICAMS:::Canonicalize1Del(context, tmp$del_str, tmp$leftmost_pos)
    prev_ret = ICAMS::Canonicalize1Del(context, ref, pos, trace) # pos is the start of the deletion

    if (is.na(prev_ret) || is.na(new_ret) || prev_ret != new_ret) {
      message("\n\nDELETION difference:")
      message(prev_ret, " vs ", new_ret, " ref = ", ref)
    }

    new_ret2 = categorize_1_indel(
      context = context,
      ins_or_del = "d",
      ins_or_del_seq = tmp$del_str,
      pos = tmp$leftmost_pos
    )
    # browser()
    if (prev_ret != new_ret2$COSMIC_83) {
      message("\n\nDELETION difference:")
      message(prev_ret, " vs ", new_ret2, " ref = ", ref)
    }
  } else if (nchar(alt) > nchar(ref)) {
    # An insertion

    tmp_short = context
    tmp_long = stringi::stri_sub_replace(
      tmp_short,
      pos,
      pos - 1, # pos was just _after_ the site of the insertion
      replacement = alt
    )

    tmp = justify_indel(tmp_long, tmp_short, pos, alt)
    new_ret = ICAMS:::Canonicalize1INS(
      context,
      tmp$del_str,
      tmp$leftmost_pos - 1
    )

    new_ret2 = categorize_1_indel(
      context = tmp_long,
      ins_or_del = "i",
      ins_or_del_seq = tmp$del_str,
      pos = tmp$leftmost_pos
    )

    prev_ret = ICAMS:::Canonicalize1INS(context, alt, pos - 1, trace) # the insertion occurs immediately after pos - 1

    if (is.na(prev_ret) || is.na(new_ret) || prev_ret != new_ret) {
      message("\n\nINSERTION difference:")
      message(prev_ret, " vs ", new_ret, " ref = ", ref)
    }
  } else {
    stop("Non-insertion / non-deletion found: ", ref, " ", alt, " ", context)
  }
  if (regress) {
    return(new_ret3$COMSMIC_83)
  } else {
    return(new_ret2)
  }
}

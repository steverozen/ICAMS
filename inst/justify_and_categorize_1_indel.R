#' @title Given a single insertion or deletion in context categorize it.
#'
#' @param context Ample surrounding
#'   sequence on each side of the insertion or deletion.
#'
#' @param orig_ref The reference allele as a single string.
#'
#' @param orig_alt The alternative allele as a single string
#'
#' @param orig_pos The position of the inserted or deleted sequence in \code{context}.
#' If this is an insertion, it is the positon after the position of the insertion of \code{alt}.
#' If this is a deletion, it is the position in context at which \code{ref} starts.
#'
#' @param explain If TRUE then generate output explaining the processing
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
justify_and_categorize_1_indel = function(
  context,
  orig_ref,
  orig_alt,
  orig_pos,
  explain_indels = FALSE,
  regress = TRUE,
  remove_common_prefix = FALSE
) {
  if (remove_common_prefix) {
    if (orig_ref == "" || orig_alt == "") {
      stop(
        "cannot remove common prefix if orig_ref or orig_alt is the empty string"
      )
    }
    stopifnot(substr(orig_ref, 1, 1) == substr(orig_alt, 1, 1))
    ref <- substr(orig_ref, 2, nchar(orig_ref))
    alt <- substr(orig_alt, 2, nchar(orig_alt))
    pos = orig_pos + 1
  } else {
    ref = orig_ref
    alt = orig_alt
    pos = orig_pos
  }

  if (nchar(alt) < nchar(ref)) {
    # We have a deletion.
    # pos is the 1-based start of the deletion in the input argumet "context".

    if (explain_indels) {
      message("\n\nExplanation =========== del of ", ref, " =====")
      message("before: ", context)
      message(
        "after:  ",
        substr(context, 1, pos - 1),
        strrep("-", nchar(ref)),
        substr(context, pos + nchar(ref), nchar(context))
      )
    }

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
    prev_ret = ICAMS::Canonicalize1Del(context, ref, pos) # pos is the start of the deletion

    if (is.na(prev_ret) || is.na(new_ret) || prev_ret != new_ret) {
      message("\n\nDELETION difference 1:")
      message(prev_ret, " vs ", new_ret, " ref = ", ref)
    }

    if (explain_indels) {
      message("\n\nExplanation =========== del of ", tmp$del_str, " =====")
      message("before: ", context)
      message(
        "after:  ",
        substr(context, 1, tmp$leftmost_pos - 1),
        strrep("-", nchar(ref)),
        substr(context, tmp$leftmost_pos + nchar(ref), nchar(context))
      )
    }

    # TATCATTTTCCATCATTCTATTCAAGCTTTTCTTCTT------TGTTACAACATTTTTGGTATTACATGACTTCTCCTA ->
    # ATCATTTTCCATCATTCTATTCAAGCTTTTCTTCTTTGTTACAACATTTTTGGTATTACATGACTTCTCCTA

    # TATCATTTTCCATCATTCTATTCAAGCTTTTC------TTCTTTGTTACAACATTTTTGGTATTACATGACTTCTCCTA ->
    # TATCATTTTCCATCATTCTATTCAAGCTTTTCTTCTTTGTTACAACATTTTTGGTATTACATGACTTCTCCTA
    # repeat is TTCTTT TTCTTT

    new_ret2 = categorize_1_justified_indel(
      context = context,
      ins_or_del = "d",
      ins_or_del_seq = tmp$del_str,
      pos = tmp$leftmost_pos
    )

    if (is.na(prev_ret) || prev_ret != new_ret2$COSMIC_83) {
      message("\n\nDELETION difference 2:")
      message("old = ", prev_ret)
      message("new = ", new_ret2$COSMIC_83, " ref = ", ref)
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

    prev_ret = ICAMS:::Canonicalize1INS(context, alt, pos - 1) # the insertion occurs immediately after pos - 1

    new_ret2 = categorize_1_justified_indel(
      context = tmp_long,
      ins_or_del = "i",
      ins_or_del_seq = tmp$del_str,
      pos = tmp$leftmost_pos
    )

    # message("\n\n\n==================\ninsertion old string = ", prev_ret)

    if (is.na(prev_ret) || is.na(new_ret) || prev_ret != new_ret) {
      message("\nINSERTION difference 1:")
      message(prev_ret, " vs ", new_ret, " alt = ", alt)
    }

    if (prev_ret != new_ret2$COSMIC_83) {
      message("\nINSERTION difference: 2")
      message("old = ", prev_ret)
      message("new = ", new_ret2$COSMIC_83, " alt = ", alt)
    }
  } else {
    stop("Non-insertion / non-deletion found: ", ref, " ", alt, " ", context)
  }
  if (regress) {
    return(new_ret2$COMSMIC_83)
  } else {
    return(new_ret2)
  }
}

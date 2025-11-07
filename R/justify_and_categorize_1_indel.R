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
#' @param explain_indels If TRUE then generate output explaining the processing
#'
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
#'
justify_and_categorize_1_indel = function(
  context,
  orig_ref,
  orig_alt,
  orig_pos,
  explain_indels = FALSE
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
    # We have a deletion.
    # pos is the 1-based start of the deletion in the input argumet "context".

    tmp_long = context
    end_pos = pos + nchar(ref) - 1
    tmp_short = stringi::stri_sub_replace(
      context,
      from = pos,
      to = end_pos,
      replacement = ""
    )

    tmp = justify_indel(tmp_long, tmp_short, pos, ref)

    explain_del = function() {
      message("\n\nExplanation =========== deletion of ", ref, " =====")
      message("Prior to justifying")
      show_indel(tmp_long, tmp_short, pos, "d")
      message("After justifying")
      show_indel(tmp_long, tmp_short, tmp$leftmost_pos, "d")
    }

    if (explain_indels) {
      explain_del()
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

    if (regress) {
      prev_ret = Canonicalize1Del(context, ref, pos) # pos is the start of the deletion

      if (is.na(prev_ret) || prev_ret != new_ret2$COSMIC_83) {
        message("\n\nDELETION difference 2:")
        message("old = ", prev_ret)
        message("new = ", new_ret2$COSMIC_83, " ref = ", ref)
        explain_del()
      }
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

    explain_ins = function() {
      message("\n\nExplanation =========== insertion of ", alt, " =====")
      message("Prior to justifying")
      show_indel(tmp_long, tmp_short, pos, "i")
      message("After justifying")
      show_indel(tmp_long, tmp_short, tmp$leftmost_pos, "i")
    }

    if (explain_indels) {
      explain_ins()
    }

    new_ret2 = categorize_1_justified_indel(
      context = tmp_long,
      ins_or_del = "i",
      ins_or_del_seq = tmp$del_str,
      pos = tmp$leftmost_pos
    )

    if (regress) {
      prev_ret = Canonicalize1INS(context, alt, pos - 1) # the insertion occurs immediately after pos - 1

      if (prev_ret != new_ret2$COSMIC_83) {
        message("\nINSERTION difference: 2")
        message("old = ", prev_ret)
        message("new = ", new_ret2$COSMIC_83, " alt = ", alt)
        explain_ins()
      }
    }
  } else {
    stop("Non-insertion / non-deletion found: ", ref, " ", alt, " ", context)
  }
  if (regress) {
    new_ret2$prev_COSMIC_83 = prev_ret
  }
  return(new_ret2)
}

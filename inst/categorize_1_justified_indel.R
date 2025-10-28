# categorize_del

library(Biostrings)

#' Given a deletion and its sequence context, categorize it
#'
#' This function is primarily for internal use, but we export it
#' to document the underlying logic.
#'
#' See \url{https://github.com/steverozen/ICAMS/blob/v3.0.9-branch/data-raw/PCAWG7_indel_classification_2021_09_03.xlsx}
#' for additional information on deletion mutation classification.
#'
#' This function first handles deletions in homopolymers, then
#' handles deletions in simple repeats with
#' longer repeat units, (e.g. \code{CACACACA}, see
#' \code{\link{FindMaxRepeatDel}}),
#' and if the deletion is not in a simple repeat,
#' looks for microhomology (see \code{\link{FindDelMH}}).
#'
#' See the code for unexported function \code{\link{CanonicalizeID}}
#' and the functions it calls for handling of insertions.
#'
#' @param context The deleted sequence plus ample surrounding
#'   sequence on each side (at least as long as \code{del.seq}).
#'
#' @param del.seq The deleted sequence in \code{context}.
#'
#' @param pos The position of \code{del.sequence} in \code{context}.
#'
#' @param trace If > 0, then generate messages tracing
#' how the computation is carried out.
#
#' @return A string that is the canonical representation
#'  of the given deletion type. Return \code{NA}
#'  and raise a warning if
#'  there is an un-normalized representation of
#'  the deletion of a repeat unit.
#'  See \code{FindDelMH} for details.
#'  (This seems to be very rare.)
#'
#' @examples
#' categorize_1_indel("xyAAAqr", "d", del.seq = "A", pos = 3) # "DEL:T:1:2"
#' categorize_1_indel("xyAAAqr", "d", del.seq = "A", pos = 4) # "DEL:T:1:2"
#' categorize_1_indel("xyAqr", "d", del.seq = "A", pos = 3)   # "DEL:T:1:0"
#'
#' @export

#  ins_or_del, previous_char, repeat_seq, repeat_count, post_char, mh, COSMIC_83

categorize_1_justified_indel <- function(
  context,
  ins_or_del,
  ins_or_del_seq,
  pos,
  trace = 0,
  verbose = 0
) {
  # Pos is the 1-based position of the first base that was deleted
  # is it 1 bp deletion?
  mh = NA
  stopifnot(pos >= 2)

  if (verbose > 0) {
    message("\n===============================")
    message("ins_or_del = ", ins_or_del)
    message("ins_or_del_seq = ", ins_or_del_seq)
    message("context = ", context)
    message("pos = ", pos)
    message("before pos context = ", substr(context, 1, pos - 1))
    message("substr(context, pos, pos) = ", substr(context, pos, pos))
    message(
      "substr(context, pos + 1, pos + 1) = ",
      substr(context, pos + 1, pos + 1)
    )
  }
  if (FALSE && nchar(ins_or_del_seq) > 3) {
    browser()
  }
  if (ins_or_del_seq == "xGGAGTGGGGCCT") {
    browser()
  }
  regex = paste0("^.{", pos - 2, "}(.)((", ins_or_del_seq, ")+)(.)(.*$)")
  mymatch = stringr::str_match(context, regex)[1, ]
  pre = mymatch[2]
  unmutated_rep_count = nchar(mymatch[3]) / nchar(ins_or_del_seq)
  if (FALSE && unmutated_rep_count > 1) {
    browser()
  }
  stopifnot(unmutated_rep_count == floor(unmutated_rep_count))
  if (ins_or_del == "i") {
    unmutated_rep_count = unmutated_rep_count - 1
    # The context arg is the sequence after the insertion. We want unmutated_rep_count to
    # reflect the repeat count prior to the mutation
  }
  post = mymatch[5]
  post_all = paste0(post, mymatch[6])
  if (verbose > 0) {
    message("regex = ", regex)
    message("unmutated_rep_count = ", unmutated_rep_count)
    message("after match")
    message("pre = ", pre)
    message("mymatch[4] (repeats) = ", mymatch[4])
    message("post = ", post)
    message("post_all = ", post_all)
  }

  if (nchar(ins_or_del_seq) == 1) {
    if (ins_or_del_seq %in% c("A", "G")) {
      pre = ICAMS::revc(post)
      ins_or_del_seq = ICAMS::revc(ins_or_del_seq)
      post = ICAMS::revc(mymatch[2]) # pre was already overwritten
    }
  } else if (unmutated_rep_count == 1) {
    if (ins_or_del == "d") {
      # Check for micrhomology
      microhomology_len = Biostrings::lcprefix(ins_or_del_seq, post_all)
      if (microhomology_len > 0) {
        mh = microhomology_len
      }
    }
  }

  retlist = list(
    ins_or_del = ins_or_del,
    pre = pre,
    ins_or_del_seq = ins_or_del_seq,
    unmutated_rep_count = unmutated_rep_count,
    post = post,
    mh = mh
  )

  retlist$COSMIC_83 = gen_COSMIC_83_string(retlist)

  return(retlist)
} # End categorize_del

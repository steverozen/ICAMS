# Take the view that there is a deletion in long_str
# at pos that creates short_str.
#
# pos is the 1-based position in long string at which to make a deletion to get short_str
# dually, pos is the 1-based position in short string before which one can
# make an insertion to get long_str.
#
# Move pos as far to left as possible so that a deletion
# at that position still results in and edit of long_str
# to short_str
#
# This can also be interpreted an inserion in short_str
# immediately in front of pos generating long_str

justify_indel = function(long_str, short_str, pos, expected_delta = NULL) {
  library(stringi)

  del_len = nchar(long_str) - nchar(short_str)

  if (!is.null(expected_delta)) {
    if (expected_delta != substr(long_str, pos, pos + del_len - 1)) {
      message(
        "expected_delta = ",
        expected_delta,
        " != substr(long_str, pos, pos + del_len - 1)) = ",
        substr(long_str, pos, pos + del_len - 1)
      )
      browser()
      stop()
    }
  }

  stopifnot(
    stri_sub_replace(
      long_str,
      from = pos,
      to = pos + del_len - 1,
      replacement = ""
    ) ==
      short_str
  )
  for (test_pos in (pos - 1):0) {
    if (
      stri_sub_replace(
        long_str,
        test_pos,
        test_pos + del_len - 1,
        replacement = ""
      ) !=
        short_str
    ) {
      break
    }
  }
  leftmost_pos = test_pos + 1
  return(
    list(
      leftmost_pos = leftmost_pos,
      del_str = stri_sub(
        long_str,
        leftmost_pos,
        leftmost_pos + del_len - 1
      )
    )
  )
}

if (FALSE) {
  justify_indel("abc", "ac", 2)
  justify_indel("abbbc", "abbc", 4)
  justify_indel("xycagcaguv", "xycaguv", 4)
  justify_indel("xycagcaguv", "xycaguv", 5)
  justify_indel("xycagcaguv", "xycaguv", 6)
  justify_indel("xycagcaguv", "xycaguv", 7) # error
  justify_indel("xycagcagcagcaguv", "xycagcagcaguv", 9)
}


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
xCanonicalize1ID <- function(context, ref, alt, pos, trace = 0) {
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
    # browser()
    new_ret = ICAMS:::Canonicalize1Del(context, tmp$del_str, tmp$leftmost_pos)
    new_ret2 = categorize_del(context, tmp$del_str, tmp$leftmost_pos)

    prev_ret = ICAMS::Canonicalize1Del(context, ref, pos, trace) # pos is the start of the deletion

    if (is.na(prev_ret) || is.na(new_ret) || prev_ret != new_ret) {
      message("\n\nDELETION difference:")
      message(prev_ret, " vs ", new_ret, " ref = ", ref)
    }
    return(prev_ret)
  } else if (nchar(alt) > nchar(ref)) {
    # An insertion

    tmp_short = context
    tmp_long = stri_sub_replace(
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

    prev_ret = ICAMS:::Canonicalize1INS(context, alt, pos - 1, trace) # the insertion occurs immediately after pos - 1

    if (is.na(prev_ret) || is.na(new_ret) || prev_ret != new_ret) {
      message("\n\nINSERTION difference:")
      message(prev_ret, " vs ", new_ret, " ref = ", ref)
    }

    return(prev_ret)
  } else {
    stop("Non-insertion / non-deletion found: ", ref, " ", alt, " ", context)
  }
}

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
#' Canonicalize1Del("xyAAAqr", del.seq = "A", pos = 3) # "DEL:T:1:2"
#' Canonicalize1Del("xyAAAqr", del.seq = "A", pos = 4) # "DEL:T:1:2"
#' Canonicalize1Del("xyAqr", del.seq = "A", pos = 3)   # "DEL:T:1:0"
#'
#' @export

categorize_del <- function(context, del.seq, pos, trace = 0) {
  # is it 1 bp deletion?
  if (nchar(del.seq) == 1) {
    browser()
    stopifnot(pos >= 2)
    regex = paste0("^.{", pos - 2, "}(.)(", del.seq, "+)([^", del.seq, "])")
    message("regex = ", regex)
    match = stringr::str_match(context, regex)[1, ]
    pre = match[2]
    rep_count = nchar(match[3])
    post = match[4]
    if (del.seq %in% c("A", "G")) {
      pre = ICAMS::revc(match[4])
      # del.seq = ICAMS::revc(del.seq)
      post = ICAMS::revc(match[2])
    }
  }

  # Is the deletion involved in a repeat?
  rep.count <- ICAMS::FindMaxRepeatDel(context, del.seq, pos)

  rep.count.string <- ifelse(rep.count >= 5, "5+", as.character(rep.count))
  deletion.size <- nchar(del.seq)
  deletion.size.string <-
    ifelse(deletion.size >= 5, "5+", as.character(deletion.size))

  # Category is "1bp deletion"
  if (deletion.size == 1) {
    if (del.seq == "G") {
      del.seq <- "C"
    }
    if (del.seq == "A") {
      del.seq <- "T"
    }
    return(paste0("DEL:", del.seq, ":1:", rep.count.string))
  }

  # Category is ">2bp deletion"
  if (rep.count > 0) {
    return(
      paste0("DEL:repeats:", deletion.size.string, ":", rep.count.string)
    )
  }

  # We have to look for microhomology
  microhomology.len <- FindDelMH(context, del.seq, pos, trace = trace)
  if (microhomology.len == -1) {
    warning(
      "Non-normalized deleted repeat ignored:",
      "\ncontext: ",
      context,
      "\ndeleted sequence: ",
      del.seq,
      "\nposition of deleted sequence: ",
      pos
    )
    return(NA)
  }
  if (microhomology.len == 0) {
    stopifnot(rep.count.string == 0)
    # Categorize and return non-repeat, non-microhomology deletion
    return(paste0("DEL:repeats:", deletion.size.string, ":0"))
  }

  microhomology.len.str <-
    ifelse(microhomology.len >= 5, "5+", as.character(microhomology.len))

  return(paste0(
    "DEL:MH:",
    deletion.size.string,
    ":",
    microhomology.len.str
  ))
}

# categorize_del

library(Biostrings)

#' Given a indel and its sequence context, categorize it
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
#' @param context The sequence surrounding the indel.
#'
#' @param ins_or_del A singgle character, with "i" denoting
#' an insertion and "d" denotine an deletion.
#'
#' @param ins_or_del_seq The the sequence that was inserted or deleted.
#'
#' @param pos For deletions, the 1-based position of the start of the
#' deleted sequence; for insertions, the position immediately to the right
#' of where the inserrtion occurs.
#'
#' @param verbose If > 0, then generate messages tracing
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
#' categorize_1_justified_indel("xyAAAqr", "d", del.seq = "A", pos = 3) # "DEL:T:1:2"
#' categorize_1_justified_indel("xyAAAqr", "d", del.seq = "A", pos = 4) # "DEL:T:1:2"
#' categorize_1_justified_indel("xyAqr", "d", del.seq = "A", pos = 3)   # "DEL:T:1:0"
#'
#' @export

#  ins_or_del, previous_char, repeat_seq, repeat_count, post_char, mh, COSMIC_83

categorize_1_justified_indel <- function(
  context,
  ins_or_del,
  ins_or_del_seq,
  pos,
  verbose = 0
) {
  # Pos is the 1-based position of the first base that was deleted
  # is it 1 bp deletion?
  mh = 0L
  koh_mh = 0L
  stopifnot(pos >= 2)
  ins_or_del_seq_len = nchar(ins_or_del_seq)

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

  regex = paste0("^.{", pos - 2, "}(.)((", ins_or_del_seq, ")+)(.)(.*$)")
  mymatch = stringr::str_match(context, regex)[1, ]
  pre = mymatch[2]

  all_repeated_seq = mymatch[3]

  # indel_str_count_in_ref is the number of times the indel string
  # appears in the reference sequence. For deletions this must be
  # >= 1, for insertions it can be 0 times.
  indel_str_count_in_ref = nchar(all_repeated_seq) / ins_or_del_seq_len
  if (ins_or_del == "i") {
    indel_str_count_in_ref = indel_str_count_in_ref - 1
    # The context arg is the sequence after the insertion. We want indel_str_count_in_ref to
    # reflect the repeat count prior to the insertion
  }

  stopifnot(indel_str_count_in_ref == floor(indel_str_count_in_ref))
  indel_str_count_in_ref = as.integer(indel_str_count_in_ref)

  post = mymatch[5]
  post_all = paste0(post, mymatch[6])
  if (verbose > 0) {
    message("regex = ", regex)
    message("indel_str_count_in_ref = ", indel_str_count_in_ref)
    message("after match")
    message("pre = ", pre)
    message("mymatch[4] (repeats) = ", mymatch[4])
    message("post = ", post)
    message("post_all = ", post_all)
  }

  if (ins_or_del_seq_len == 1) {
    R = indel_str_count_in_ref
    U = 1L
    if (ins_or_del_seq %in% c("A", "G")) {
      pre = ICAMS::revc(post)
      ins_or_del_seq = ICAMS::revc(ins_or_del_seq)
      post = ICAMS::revc(mymatch[2]) # pre was already overwritten
    }
  } else {
    # For the Koh et al. 2025 classification we need to see if there
    # are repeats within ins_or_del_seq.  See Fig 2 A from this paper.

    # newpattern will match the shortest prefix, p, of a string, x,
    # such that x = 0 or more repeats6 of p. One note: if we would
    # say "there are no repeats in ins_or_del_seq", then p is
    # exactly ins_or_del_seq.

    newpattern = "^(.+?)\\1*$"
    newmatch = stringr::str_match(ins_or_del_seq, newpattern)
    shortest_prefix = newmatch[1, 2]

    # U if from the nomenclature in Koh et al, Fig 2a.
    U = nchar(shortest_prefix)

    # Is shortest_prefix repeated in post_all?
    R_match_pattern = paste0("^(?:", shortest_prefix, ")+")
    R_match = stringr::str_match(
      paste0(all_repeated_seq, post_all),
      R_match_pattern
    )

    # We are using nomenclature from Koh et al. again here.
    R = (nchar(R_match[1, 1]) / U)
    if (is.na(R)) {
      browser() # This is a programming error
    }

    R_outside_ins_or_del_seq = R - (ins_or_del_seq_len / U)
    if (ins_or_del == "i") {
      R = R_outside_ins_or_del_seq
    } else {
      if (R_outside_ins_or_del_seq == 0) {
        # browser()
        # I think we need to set U to ins_or_del_seq_len and set R to 1 See lines
      }
    }

    stopifnot(R == floor(R))
    R = as.integer(R)

    if (ins_or_del == "d") {
      if (indel_str_count_in_ref == 1) {
        # Check for micrhomology based on the ins_or_del_seq alone
        mh = Biostrings::lcprefix(ins_or_del_seq, post_all)
        if (length(R) == 0) {
          browser() # This is programming error
        }
        if (R == 1) {
          koh_mh = mh
        }
      }
    } else {
      # Insertion
      if (indel_str_count_in_ref == 0) {
        mh = Biostrings::lcprefix(ins_or_del_seq, post_all)
        if (R == 0) {
          # shortest_prefix can be shorted than ins_or_del_seq, for example
          # in the insertion ATC|GG|TC where GG is inserted.
          # Then shortest_prefix is G but ins_or_del_seq is GG.

          koh_mh = mh
        }
      }
    }
  } # end else (i.e. ins_or_del_seq_len > 1)

  retlist = list(
    ins_or_del = ins_or_del,
    pre = pre,
    ins_or_del_seq = ins_or_del_seq,
    indel_str_count_in_ref = indel_str_count_in_ref,
    post = post,
    mh = mh,
    R = R,
    U = U,
    koh_mh = koh_mh
  )

  retlist$COSMIC_83 = ICAMS:::gen_COSMIC_83_string(retlist)
  retlist$Koh_89 = gen_Koh_89_string(retlist)
  retlist$Koh_476 = gen_Koh_476_string(retlist)

  return(retlist)
} # End categorize_del

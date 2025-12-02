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
#' @param context The sequence surrounding the indel PRIOR to the insertion
#' or deletion.
#'
#' @param ins_or_del A singgle character, with "i" denoting
#' an insertion and "d" denotine an deletion.
#'
#' @param ins_or_del_seq The the sequence that was inserted or deleted.
#'
#' @param pos For deletions, the 1-based position of the start of the
#' deleted sequence; for insertions, the position immediately to the right
#' of where the insertion occurs.
#'
#' @return A string that is the canonical representation
#'  of the given deletion type. Return \code{NA}
#'  and raise a warning if
#'  there is an un-normalized representation of
#'  the deletion of a repeat unit.
#'  See \code{FindDelMH} for details.
#'  (This seems to be very rare.)
#'
#' @examples
#' categorize_1_justified_indel("GGAAAGG", "d", ins_or_del_seq = "A", pos = 3) # "DEL:T:1:2"
#' categorize_1_justified_indel("GGAAAGG", "d", ins_or_del_seq = "A", pos = 4) # "DEL:T:1:2"
#' categorize_1_justified_indel("TTATT", "d", ins_or_del_seq = "A", pos = 3)   # "DEL:T:1:0"
#'

categorize_1_justified_indel <- function(
  context,
  ins_or_del,
  ins_or_del_seq,
  pos
) {
  if (TRUE) {
    if (ins_or_del == "i") {
      start_of_slice3 = pos
    } else {
      start_of_slice3 = pos + nchar(ins_or_del_seq)
    }
    slice3 = stringi::stri_sub(context, from = start_of_slice3)
    if (ins_or_del == "d") {
      if (
        paste0(ins_or_del_seq, slice3) != stringi::stri_sub(context, from = pos)
      ) {
        browser()
      }
    }
    koh_extra = seg_simple(
      ins_or_del = ins_or_del,
      string = ins_or_del_seq,
      context = slice3
    )
  } else {
    koh_extra = list()
  }

  mh = 0L
  koh_mh = 0L
  if (pos < 2) {
    message(
      "position of the insertion or deletion of ",
      ins_or_del_seq,
      " is ",
      pos
    )
    message("This should be >= 2")
    err_ret = indel_all_na_return()
    return(err_ret)
  }
  ins_or_del_seq_len = nchar(ins_or_del_seq)

  if (grepl("N", ins_or_del_seq)) {
    message(
      "N found in indel sequence: ",
      ins_or_del_seq,
      " context: ",
      context
    )
  }

  regex = paste0("^.{", pos - 2, "}(.)((", ins_or_del_seq, ")+)(.)(.*$)")

  if (ins_or_del == "i") {
    x_context = stringi::stri_sub_replace(
      context,
      pos,
      pos - 1, # pos was just _after_ the site of the insertion
      replacement = ins_or_del_seq
    )
  } else {
    x_context = context
  }

  mymatch = stringr::str_match(x_context, regex)[1, ]
  pre = mymatch[2]

  all_repeated_seq = mymatch[3]

  # indel_str_count_in_ref is the number of times the indel string
  # appears in the reference sequence prior to the mutations.
  # Therefore, for deletions this must be >= 1, for insertions it can be 0 times.
  indel_str_count_in_ref = nchar(all_repeated_seq) / ins_or_del_seq_len
  if (ins_or_del == "i") {
    indel_str_count_in_ref = indel_str_count_in_ref - 1
    # x_context was the sequence after the insertion. We want
    # indel_str_count_in_ref to reflect the repeat count prior to the insertion
  }

  stopifnot(indel_str_count_in_ref == floor(indel_str_count_in_ref))

  post = mymatch[5]
  post_all = paste0(post, mymatch[6])

  if (ins_or_del_seq_len == 1) {
    R = indel_str_count_in_ref
    U = 1L
    if (ins_or_del_seq %in% c("A", "G")) {
      pre = ICAMS::revc(post)
      ins_or_del_seq = ICAMS::revc(ins_or_del_seq)
      post = ICAMS::revc(mymatch[2]) # pre was already overwritten
    }
    U_seq = ins_or_del_seq
    U_seq_count_in_indel_seq = 1
    R_outside_ins_or_del_seq = nchar(post_all)
  } else {
    # For the Koh et al. 2025 classification we need to see if there
    # are repeats within ins_or_del_seq.  See Fig 2 A from this paper.

    # newpattern will match the shortest prefix, p, of a string, x,
    # such that x = 0 or more repeats6 of p. One note: if we would
    # say "there are no repeats in ins_or_del_seq", then p is
    # exactly ins_or_del_seq.

    newpattern = "^(.+?)\\1*$"
    newmatch = stringr::str_match(ins_or_del_seq, newpattern)
    U_seq = newmatch[1, 2]

    # U if from the nomenclature in Koh et al, Fig 2a.
    U = nchar(U_seq)
    U_seq_count_in_indel_seq = ins_or_del_seq_len / U
    stopifnot(U_seq_count_in_indel_seq == floor(U_seq_count_in_indel_seq))

    # Is U_seq repeated in post_all?
    R_match_pattern = paste0("^(?:", U_seq, ")+")
    R_match = stringr::str_match(
      paste0(all_repeated_seq, post_all),
      R_match_pattern
    )

    # We are using nomenclature from Koh et al. again here.
    R = (nchar(R_match[1, 1]) / U)
    if (is.na(R)) {
      browser() # This is a programming error
    }

    R_outside_ins_or_del_seq = R - U_seq_count_in_indel_seq
    if (ins_or_del == "i") {
      R = R_outside_ins_or_del_seq
    } else {
      if (R_outside_ins_or_del_seq == 0) {
        # browser()
        # Do we need to set U to ins_or_del_seq_len and set R to 1?
      }
    }

    stopifnot(R == floor(R))

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
          # mh can be shorter than ins_or_del_seq, for example
          # in the insertion ATC|GG|GTC where GG is inserted.
          # Then mh is G but ins_or_del_seq is GG.

          koh_mh = mh
        }
      }
    }
  } # end else (i.e. ins_or_del_seq_len > 1)

  retlist = c(
    list(
      ins_or_del = ins_or_del,
      pre = pre,
      ins_or_del_seq = ins_or_del_seq,
      post = post,
      L = as.integer(ins_or_del_seq_len),
      U_seq = U_seq,
      U = as.integer(U),
      U_seq_count_in_indel_seq = as.integer(U_seq_count_in_indel_seq),
      indel_str_count_in_ref = as.integer(indel_str_count_in_ref),
      R = as.integer(R),
      R_outside_ins_or_del_seq = as.integer(R_outside_ins_or_del_seq),
      mh = as.integer(mh),
      koh_mh = as.integer(koh_mh)
    ),
    koh_extra
  )

  retlist$COSMIC_83 = gen_COSMIC_83_string(retlist)
  retlist$Koh_89 = gen_Koh_89_string(retlist)
  retlist$Koh_476 = gen_Koh_476_string(retlist)

  return(retlist)
} # End categorize_del


indel_all_na_return = function(info_string = "Unable_to_categorize") {
  retlist = list(
    ins_or_del = NA,
    pre = NA,
    ins_or_del_seq = NA,
    post = NA,
    L = NA,
    U_seq = NA,
    U = NA,
    U_seq_count_in_indel_seq = NA,
    indel_str_count_in_ref = NA,
    R = NA,
    R_outside_ins_or_del_seq = NA,
    mh = NA,
    koh_mh = NA,
    COSMIC_83 = info_string,
    Koh_89 = info_string,
    Koh_476 = info_string
  )
  return(retlist)
}

xgen_Koh_476_string = function(arglist) {
  if (arglist$ins_or_del == "d") {
    INS_OR_DEL = "Del"
  } else {
    stopifnot(arglist$ins_or_del == "i")
    INS_OR_DEL = "Ins"
  }

  rep_count = arglist$R

  ins_or_del_seq = arglist$ins_or_del_seq
  ins_or_del_len = nchar(ins_or_del_seq)

  if (ins_or_del_len == 1) {
    # Lines 4 through 183 (insertions) and 244 through 405 (deletions) of Koh et al. sup table 7
    repcountstr = ifelse(rep_count >= 9, "9+", as.character(rep_count))
    return(paste0(
      arglist$pre,
      "[",
      INS_OR_DEL,
      "(",
      ins_or_del_seq,
      "):R",
      repcountstr,
      "]",
      arglist$post
    ))
  }

  microhom_len = arglist$mh_koh

  if (!is.na(microhom_len && microhom_len > 0)) {
    if (INS_OR_DEL == "Ins") {
      # Insertion with microhomology
      # Lines 184 and 185
      if (rep_count != 0) {
        browser()
      }

      return(
        paste0(
          INS_OR_DEL,
          ifelse(ins_or_del_seq <= 4, "(2,4):M", "(5,):M")
        )
      )
    } else {
      # Deletion with microhomology
      # Lines 454 through 474
      if (rep_count != 1) {
        browser()
      }

      del_mh_str = ifelse(microhom_len >= 6, "(6,)", microhom_len)
      del_len_str = ifelse(ins_or_del_len >= 7, "(7,)", ins_or_del_len)
      return(paste0(INS_OR_DEL, del_len_str, ":M", del_mh_str))
    }
  }

  # Remaining classes of insertions and deletions;
  #
  # For deletion, L is the total length of the (now deleted)
  # sequence, U is the length of the repeat unit (and there might be
  # multiple repeat units within the deleted seqeucne), and R is
  # the numbuer of repeat units prior to the deletion.

  # For insertions, L is the total length of the inserted sequence
  # and R is the number of repeat units prior to the
  # insertion. We don't use U in the description, and figure 2
  # in Koh et al doesn't show an example of cases where the
  # inerted sequence has multiple repeat units, but we
  # assume that R is defiened analogously to the case for
  # insertions, so R is the number of repeat units

  # There is an edge case:
  #
  # |ACAC|ACGTG where the intital sequence is either an insertion or deletion.
  #
  # If this is a a deletion this could be classified as microhomology Del4:M2.
  # But could also be classified as Del4:U2:R3 (Koh-476) or Del(2,8):U(1,2)R(2:4) (Koh-89)
  #
  # If this is an insertion it could be Ins4:U2:R1 or Ins(2,4):M.  Since
  # all Ins4:U2:R1 could be considered Ins(2,4):M, we assume that
  # Ins4:U2:R1 takes precedence, and the analgous reasoning applies to
  # deletions.

  # We as assume that for e.g. |ABABA|BABAB we consider U = 5, L = 5, and for insertion R = 0

  tt = function(s) {
    pattern <- "^(.+?)\\1*$"
    r1 <- stringr::str_match(s, pattern)
    shortest_prefix = r1[1, 2]
    return(shortest_prefix)
  }

  rep_count_string = ifelse(rep_count >= 5, "5+", rep_count)
  rm(rep_count)
  if (nchar(arglist$ins_or_del_seq) == 1) {
    return(paste0(INS_OR_DEL, arglist$ins_or_del_seq, ":1:", rep_count_string))
  }

  size_string = nchar(arglist$ins_or_del_seq)
  if (size_string >= 5) {
    size_string = "5+"
  } else {
    size_string = as.character(size_string)
  }

  if (arglist$ins_or_del == "d" && !is.na(microhom_len) && microhom_len > 0) {
    mh_string = ifelse(microhom_len >= 5, "5+", microhom_len)
    return(paste0("DEL:MH:", size_string, ":", mh_string))
  }

  return(paste0(
    INS_OR_DEL,
    "repeats:",
    size_string,
    ":",
    rep_count_string
  ))
}

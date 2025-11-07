gen_Koh_476_string = function(arglist) {
  if (arglist$ins_or_del == "d") {
    INS_OR_DEL = "Del"
  } else {
    stopifnot(arglist$ins_or_del == "i")
    INS_OR_DEL = "Ins"
  }

  R = arglist$R

  ins_or_del_seq = arglist$ins_or_del_seq
  L = nchar(ins_or_del_seq)

  if (L == 1) {
    if (!ins_or_del_seq %in% c("A", "C", "G", "T")) {
      return(paste0("Cannot categorize indel of ", ins_or_del_seq))
    }
    if (!arglist$pre %in% c("A", "C", "G", "T")) {
      return(paste0("Cannot categorize indel preceded by  ", arglist$pre))
    }
    if (!arglist$post %in% c("A", "C", "G", "T")) {
      return(paste0("Cannot categorize indel followed by  ", arglist$post))
    }
    # Lines 4 through 183 (insertions) and 244 through 405 (deletions) of Koh et al. sup table 7
    R_str = ifelse(R >= 9, "9+", as.character(R))
    return(paste0(
      arglist$pre,
      "[",
      INS_OR_DEL,
      "(",
      ins_or_del_seq,
      "):R",
      R_str,
      "]",
      arglist$post
    ))
  }

  microhom_len = arglist$koh_mh

  if (microhom_len > 0) {
    if (INS_OR_DEL == "Ins") {
      # Insertion with microhomology
      # Lines 184 and 185
      if (R != 0) {
        browser() # This should be an error (?)
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
      if (R != 1) {
        browser() # This should be an error (?)
      }

      del_mh_str = ifelse(microhom_len >= 6, "(6,)", microhom_len)
      del_len_str = ifelse(L >= 7, "(7,)", L)
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

  U = arglist$U

  if (INS_OR_DEL == "Ins") {
    L_str = ifelse(L >= 5, "(5,)", L)
    R_str = ifelse(R >= 5, "(5,)", R)
    if (R == 0) {
      if (L >= 5) {
        return("Ins(5,):R0")
      } else {
        U_str = as.character(U)
      }
    } else {
      U_str = ifelse(U >= 3, "(3,)", U)
    }
  } else {
    L_str = ifelse(L >= 6, "(6,)", L)
    U_str = ifelse(U >= 4, "(4,)", U)
    R_str = ifelse(R >= 7, "(7+)", R)

    if (L >= 6) {
      if (U == 1) {
        return("Del(6,):U1:R(7,)")
      }
      if (U == 2) {
        return("Del(6,):U2:R(4,)")
      }
      if (U == 3) {
        return("Del(6,):U3:R(3,)")
      }
      if (U >= 4) {
        return("Del(6,):U(4,):R(2,)")
      }
    }
    if (L == 5 && U == 1) {
      return("Del5:U1:R(6,9)")
    }
  }

  # Del(6,):U2:R(4,9)

  paste0(INS_OR_DEL, L_str, ":U", U_str, ":R", R_str)
}

gen_Koh_476_string = function(arglist) {
  open_interval_format = FALSE
  fiveplus_str = ifelse(open_interval_format, "(5,)", "(5,9)")

  if (arglist$ins_or_del == "d") {
    INS_OR_DEL = "Del"
  } else {
    stopifnot(arglist$ins_or_del == "i")
    INS_OR_DEL = "Ins"
  }

  R = arglist$R

  ins_or_del_seq = arglist$ins_or_del_seq
  L = arglist$L
  U = arglist$U

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
    R_str = ifelse(
      R >= 13,
      ifelse(open_interval_format, "13+", "13"),
      as.character(R)
    )
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
  # arglist$indel_str_count_in_ref
  if (
    (INS_OR_DEL == "Del" && arglist$R_outside_ins_or_del_seq == 0 && U != 1) || # not sure about U != 1 I think Koh does not have this...
      ## (INS_OR_DEL == "Del" && arglist$indel_str_count_in_ref == 1) ||
      (INS_OR_DEL == "Ins" && arglist$indel_str_count_in_ref == 0 && U != 1)
  ) {
    microhom_len = arglist$mh
    # message("mh")

    if (INS_OR_DEL == "Del") {
      R = 1
    } else {
      R = 0
    }
  } else {
    microhom_len = arglist$koh_mh
    # message("hmh")
  }

  #  Hyptothesis: If U is 1 koh doesn't call Microhomology, s

  if (microhom_len > 0) {
    # was micohom_len
    if (INS_OR_DEL == "Ins") {
      # Insertion with microhomology
      # Lines 184 and 185
      if (R != 0) {
        # browser() # This should be an error (?)
      }
      if (L >= 5) {
        return("Ins(5,):M")
      } else {
        return("Ins(2,4):M")
      }
    } else {
      # Deletion with microhomology
      # Lines 454 through 474
      if (R != 1) {
        # browser() # This should be an error (?)
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

  if (INS_OR_DEL == "Ins") {
    L_str = ifelse(L >= 5, "(5,)", L)
    R_str = ifelse(R >= 5, fiveplus_str, R)
    if (R == 0) {
      if (L >= 5) {
        return("Ins(5,):R0")
      } else {
        return(paste0("Ins", L, ":U", U, ":R0"))
      }
    } else {
      if (L >= 5) {
        U_str = ifelse(U >= 3, "(3,)", U)
        return(paste0("Ins(5,):U", U_str, ":R", R_str))
      } else {
        return(paste0("Ins", L, ":U", U, ":R", R_str))
      }
    }
  } else {
    stopifnot(INS_OR_DEL == "Del")

    if (R == 1) {
      L_str = ifelse(L >= 10, "(10,)", L)
      if (U == 1) {
        return(paste0("Del", L_str, ":U1:R1"))
      } else {
        return(paste0("Del", L_str, ":U(2,):R1"))
      }
    }

    if (L %in% 2:4) {
      R_str = ifelse(R >= 5, fiveplus_str, R)
      return(paste0("Del", L, ":U", U, ":R", R_str))
    }

    if (L == 5) {
      if (U == 1) {
        if (R < 5) {
          return(paste0("*Del5:U1:R", R))
        }
        return("Del5:U1:R(5,9)")
      } else {
        R_str = ifelse(R >= 5, fiveplus_str, R)
        return(paste0("Del", L, ":U", U, ":R", R_str))
      }
    }

    stopifnot(L >= 6)
    U_str = ifelse(U >= 5, "(5,)", U)
    R_str = ifelse(R >= 5, ifelse(open_interval_format, "(5,)", "(5,9)"), R)

    if (U == 1) {
      if (open_interval_format) {
        return("Del(6,):U1:R(7,)")
      } else {
        return("Del(6,):U1:R(7,9)")
      }
    }
    if (U == 2) {
      if (open_interval_format) {
        return("Del(6,):U2:R(4,)")
      } else {
        return("Del(6,):U2:R(4,9)")
      }
    }
    if (U == 3) {
      if (open_interval_format) {
        return("Del(6,):U3:R(3,)")
      } else {
        return("Del(6,):U3:R(3,9)")
      }
    }
    if (U >= 4) {
      if (open_interval_format) {
        return("Del(6,):U(4,):R(2,)")
      } else {
        return(("Del(6,):U(4,):R(2,9)"))
      }
    }
  }
  stop("Should not get here: programming error")
}

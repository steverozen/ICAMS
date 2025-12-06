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
    R_str = ifelse(
      R >= 99, # Mo, adjust accoring to what you find in the actual data
      "(99,)",
      R
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

  if (arglist$spacer_length > 0 && arglist$prime3_reps > 0) {
    return(gen_koh476_mh_str(arglist))
  }

  if (
    arglist$unit_length == 1 && arglist$prime3_reps == 0 && INS_OR_DEL == "Del"
  ) {
    L_str = ifelse(L >= 10, "(10,)", L)
    return(as.character(glue::glue("Del{L_str}:U1:R1")))
  }

  if (INS_OR_DEL == "Ins") {
    # browser()
    L_str = ifelse(L >= 5, "(5,)", L)
    R_str = ifelse(arglist$R >= 5, fiveplus_str, arglist$R)
    if (arglist$R == 0) {
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

    if (arglist$spacer_length == 0 && arglist$prime3_reps == 0) {
      U_str = ifelse(arglist$unit_length >= 2, "(2,)", "1")
      if (L >= 10) {
        return(as.character(glue::glue("Del(10,):U{U_str}:R1")))
      } else {
        return(as.character(glue::glue("Del{L}:U{U_str}:R1")))
      }
    }

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

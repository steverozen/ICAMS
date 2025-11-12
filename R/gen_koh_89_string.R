gen_Koh_89_string = function(arglist) {
  if (arglist$ins_or_del == "d") {
    INS_OR_DEL = "Del"
  } else {
    stopifnot(arglist$ins_or_del == "i")
    INS_OR_DEL = "Ins"
  }

  R = arglist$R
  ins_or_del_seq = arglist$ins_or_del_seq
  L = nchar(ins_or_del_seq)
  pre = arglist$pre
  post = arglist$post

  if (L == 1) {
    if (!ins_or_del_seq %in% c("A", "C", "G", "T")) {
      return(paste0("Cannot categorize indel of ", ins_or_del_seq))
    }
    if (!pre %in% c("A", "C", "G", "T")) {
      return(paste0("Cannot categorize indel preceded by  ", pre))
    }
    if (!post %in% c("A", "C", "G", "T")) {
      return(paste0("Cannot categorize indel followed by  ", post))
    }

    if (INS_OR_DEL == "Ins") {
      # Lines 4 through 35 of Koh et al. sup table 6
      if (ins_or_del_seq == "C") {
        if (pre == "A" && R == 0 && post %in% c("A", "T")) {
          return(paste0("A[Ins(C):R0]", post))
        } else {
          if (R %in% 0:3) {
            return("Ins(C):R(0,3)")
          } else if (R %in% 4:6) {
            return("Ins(C):R(4,6)")
          } else {
            return("Ins(C):R(7,)")
          }
        }
      } else if (ins_or_del_seq == "T") {
        if (R %in% 0:4) {
          return(paste0(pre, "[Ins(T):R(0,4)]", post))
        }
        if (R %in% 5:7) {
          return(paste0(pre, "[Ins(T):R(5,7)]", post))
        }
        return(paste0(pre, "[Ins(T):R(8,)]", post))
      } else {
        browser() # an error
      }
    } else if (INS_OR_DEL == "Del") {
      if (ins_or_del_seq == "C") {
        if (R >= 6) {
          return("Del(C):R(6,9)") ## Huh? no opening bracket here, but backet for R(1,5)
        } else if (post == "G") {
          return("[Del(C):R(1,5)]G")
        }
        R_str = ifelse(R >= 4, "(4,5)", R)
        return(paste0("[Del(C):R", R_str, "]", post))
      } else if (ins_or_del_seq == "T") {
        if (R %in% 1:4) {
          R_str = "1,4"
        } else if (R %in% 5:7) {
          R_str = "5,7"
        } else {
          R_str = "8,"
        }
        return(paste0(pre, "[Del(T):R(", R_str, ")]", post))
      } else {
        browser() # A programming error
      }
    }
  }

  # L > 1 ################################

  if (INS_OR_DEL == "Ins") {
    # browser()
    indel_seq_count_in_ref = arglist$indel_str_count_in_ref

    if (FALSE) {
      # This gives us 13 deltas where we call Ins(2,4):R0
      # Without this we get 9 deltas where Koh calls Ins(5,):0
      if (arglist$mh > 0) {
        if (L <= 4) {
          return("Ins(2,4):R0")
        } else {
          return("Ins(5,):R0")
        }
      }
    }

    testR = R
    # It looks like, for insertions, the Koh paper takes R to be the indel_seq_count_in_ref
    if (testR <= 1) {
      if (L <= 4) {
        return(paste0("Ins(2,4):R", testR))
      } else {
        return(paste0("Ins(5,):R", testR))
      }
    } else {
      if (testR <= 4) {
        return("Ins(2,):R(2,4)")
      } else {
        return("Ins(2,):R(5,)")
      }
    }
  }

  ### L > 1 && Del ######################

  stopifnot(INS_OR_DEL == "Del")

  whole_del_microhom_len = arglist$mh

  if (whole_del_microhom_len > 0 && FALSE) {
    # Deletion with microhomology
    if (R != 1) {
      browser() # This would be a programming error
    }
    if (L <= 5) {
      mh_str = ifelse(
        whole_del_microhom_len <= 3,
        whole_del_microhom_len,
        "(4,)"
      )
      return(paste0("Del(2,5):M", mh_str))
    }
    mh_str = ifelse(whole_del_microhom_len <= 3, whole_del_microhom_len, "(4,)")
    return(paste0("Del(6,):M", mh_str))
  }

  koh_microhom_len = arglist$koh_mh

  if (koh_microhom_len > 0) {
    # Deletion with microhomology
    if (R != 1) {
      browser()
    }
    if (L %in% 2:5 && koh_microhom_len == 1) {
      return("Del(2,5):M1")
    } else if (L %in% 3:5 && koh_microhom_len == 2) {
      return("Del(3,5):M2")
    } else if (L %in% 4:5) {
      return("Del(4,5):M(3,4)")
    } else {
      stopifnot(L >= 6)
      if (koh_microhom_len >= 4) {
        return("Del(6,):M(4,)")
      } else {
        return(paste0("Del(6,):M", koh_microhom_len))
      }
    }
  } # koh_microhom_len > 0

  if (R == 1) {
    if (L <= 4) {
      return("Del(2,4):R1")
    }
    return("Del(5,):R1")
  }
  U = arglist$U

  if (U <= 2) {
    if (R <= 4) {
      return("Del(2,):U(1,2):R(2,4)")
    }
    return("Del(2,):U(1,2):R(5,)")
  }

  if (R == 2) {
    return("Del(3,):U(3,):R2")
  }

  "Del(3,):U(3,):R(3,)"
}

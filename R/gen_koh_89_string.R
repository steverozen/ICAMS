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
        if (pre == "A" && R == 0) {
          if (post == "A") {
            return("A[Ins(C):R0]A")
          }
          if (post == "T") {
            return("A[Ins(C):R0]T")
          }
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
        if (R %in% 5:6) {
          return(paste0(pre, "[Ins(T):R(5,6)]", post))
        }
        if (R %in% 7:8) {
          return(paste0(pre, "[Ins(T):R(7,8)]", post))
        }
        return(paste0(pre, "[Ins(T):R(9,)]", post))
      } else {
        browser() # an error
      }
    } else if (INS_OR_DEL == "Del") {
      if (ins_or_del_seq == "C") {
        if (R >= 6) {
          return("Del(C):R(6,9)")
        }
        if (post == "G") {
          return("Del(C):R(1,5)]G")
        }
        R_str = ifelse(R >= 4, "(4,5)", R)
        return(paste0(pre, "[Del(C):R", R_str, "]", post))
      } else if (ins_or_del_seq == "T") {
        if (R %in% 1:4) {
          R_str = "1,4"
        }
        if (R %in% 5:7) {
          R_str = "5,7"
        }
        R_str = "(8,)"
        return(paste0(pre, "[Del(T):R", R_str, "]", post))
      } else {
        browser() # A programming error
      }
    }
  }

  # L > 1

  if (INS_OR_DEL == "Ins") {
    if (R <= 1) {
      if (L <= 4) {
        return(paste0("Ins(2,4):R", R))
      }
      return(paste0("Ins(5,);R", R))
    }
    if (R < 4) {
      return("Ins(2,):R(2,4")
    }
    return("Ins(2,):R(5,)")
  }
  stopifnot(INS_OR_DEL == "Del")

  microhom_len = arglist$koh_mh

  if (microhom_len > 0) {
    # Deletion with microhomology
    if (R != 1) {
      browser()
    }
    if (L <= 4) {
      mh_str = ifelse(microhom_len <= 2, microhom_len, "(3,4)")
      return(paste0("del(2,", microhom_len, "):M", mh_str))
    }
    mh_str = ifelse(microhom_len <= 3, microhom_len, "(4,))")
    return(paste0("Del(6,):M", mh_str))
  }

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
    return("Del(3,):U(3):R2")
  }

  "Del(3,):U(3,):R(3,)"
}

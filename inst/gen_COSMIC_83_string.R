gen_COSMIC_83_string = function(arglist) {
  if (arglist$ins_or_del == "d") {
    INS_OR_DEL = "DEL:"
    rep_count = arglist$rep_count - 1
  } else {
    stopifnot(arglist$ins_or_del == "i")
    INS_OR_DEL = "INS:"
    rep_count = arglist$rep_count
  }
  rep_count_string = ifelse(rep_count >= 5, "5+", rep_count)
  rm(rep_count)
  if (nchar(arglist$ins_or_del_seq) == 1) {
    return(paste0(INS_OR_DEL, arglist$ins_or_del_seq, ":1:", rep_count_string))
  } else {
    if (arglist$ins_or_del == "d") {
      size_string = as.character(nchar(arglist$ins_or_del_seq))
      if (size_string >= 5) {
        size_string = "5+"
      }
      if (!is.na(arglist$mh)) {
        mh_string = ifelse(arglist$mh >= 5, "5+", arglist$mh)
        return(paste0("DEL:MH:", arglist$mh, size_string))
      } else {
        return(paste0(
          INS_OR_DEL,
          "repeats:",
          size_string,
          ":",
          rep_count_string - 1
        ))
      }
    }
  }
}

gen_Koh_476_string = function(arglist) {
  if (arglist$ins_or_del == "d") {
    INS_OR_DEL = "Del"
    rep_count = arglist$unmutated_rep_count - 1
  } else {
    stopifnot(arglist$ins_or_del == "i")
    INS_OR_DEL = "Ins"
    rep_count = arglist$unmutated_rep_count
  }

  ins_or_del_str = args$ins_or_del_str

  if (nchar(ins_or_del_str) == 1) {
    repcountstr = ifelse(rep_count >= 9, "9+", as.character(rep_count))
    return(paste0(
      arglist$pre,
      "[",
      INS_OR_DEL,
      "(",
      ins_or_del_str,
      ")R]",
      repcountstr,
      arglist$post
    ))
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

  if (arglist$ins_or_del == "d" && !is.na(arglist$mh)) {
    mh_string = ifelse(arglist$mh >= 5, "5+", arglist$mh)
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

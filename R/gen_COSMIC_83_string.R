gen_COSMIC_83_string = function(arglist) {
  if (arglist$ins_or_del == "d") {
    INS_OR_DEL = "DEL:"
    rep_count = arglist$indel_str_count_in_ref - 1
  } else {
    stopifnot(arglist$ins_or_del == "i")
    INS_OR_DEL = "INS:"
    rep_count = arglist$indel_str_count_in_ref
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

  if (arglist$ins_or_del == "d" && !is.na(arglist$mh) && arglist$mh > 0) {
    mh_string = ifelse(arglist$mh >= 5, "5+", arglist$mh)
    return(paste0("DEL:MH:", size_string, ":", mh_string))
  }

  paste0(
    INS_OR_DEL,
    "repeats:",
    size_string,
    ":",
    rep_count_string
  )
}

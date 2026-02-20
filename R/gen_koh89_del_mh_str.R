gen_koh89_del_mh_str = function(arglist) {
  stopifnot(arglist$spacer_length > 0)
  # This will be a microhomology category as
  # determined by Koh et al.'s "segment" analysis.
  stopifnot(arglist$prime_rep > 0)
  seg_microhom_len = min(
    arglist$unit_length + nchar(arglist$internal_rep),
    nchar(arglist$prime3_rep)
  )

  if (arglist$L >= 6) {
    L_str = "6,"
    if (seg_microhom_len >= 4) {
      microhom_len_str = "(4,)"
    } else {
      microhom_len_str = as.character(seg_microhom_len)
    }
  } else {
    L_str = paste0(min(4, seg_microhom_len + 1), ",5")
    if (seg_microhom_len >= 3 && seg_microhom_len <= 4) {
      microhom_len_str = "(3,4)"
    } else if (seg_microhom_len >= 6) {
      microhom_len_str = "(6,)"
    } else {
      microhom_len_str = as.character(seg_microhom_len)
    }
  }
  return(paste0("Del(", L_str, "):M", microhom_len_str))
}

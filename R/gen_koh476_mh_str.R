gen_koh476_mh_str = function(arglist) {
  stopifnot(arglist$spacer_length > 0)
  # This will be a microhomology category as
  # determined by Koh et al.'s "segment" analysis.
  stopifnot(arglist$prime_rep > 0)
  seg_microhom_len = min(
    arglist$unit_length + nchar(arglist$internal_rep),
    nchar(arglist$prime3_rep)
  )

  is_ins = arglist$ins_or_del == "i"

  if (is_ins) {
    if (arglist$L >= 5) {
      return("Ins(5,):M")
    } else {
      return("Ins(2,4):M")
    }
  }

  L_str = if (arglist$L >= 7) "(7,)" else arglist$L
  M_str = if (seg_microhom_len >= 6) "(6,)" else seg_microhom_len

  return(paste0("Del", L_str, ":M", M_str))
}

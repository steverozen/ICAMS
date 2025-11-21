explain_indel_justification = function(
  long_str,
  short_str,
  pos_before,
  pos_after,
  is_ins
) {
  len_diff = nchar(long_str) - nchar(short_str)

  message(
    "\n\nJustification explanation: ",
    ifelse(is_ins, "insertion", "deletion"),
    " of ",
    stringi::stri_sub(
      long_str,
      pos_before,
      pos_before + len_diff - 1
    ),
    " ========="
  )
  i_or_d = ifelse(is_ins, "i", "d")
  message("Prior to justifying")
  show_indel(long_str, short_str, pos_before, i_or_d)
  message("After justifying")
  show_indel(long_str, short_str, pos_after, i_or_d)
  message(
    "After justifying the ",
    ifelse(is_ins, "inserted", "deleted"),
    " sequence is ",
    stringi::stri_sub(
      long_str,
      pos_after,
      pos_after + len_diff - 1
    )
  )
}

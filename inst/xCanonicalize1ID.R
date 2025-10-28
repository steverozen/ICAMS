# temporary file to be called from ICAMS

xCanonicalize1ID = function(
  context,
  orig_ref,
  orig_alt,
  orig_pos,
  explain = FALSE,
  regress = TRUE
) {
  justify_and_categorize_1_indel(
    context,
    orig_ref,
    orig_alt,
    orig_pos,
    explain = FALSE,
    regress = TRUE,
    remove_common_prefix = FALSE
  )
}

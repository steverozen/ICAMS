prob_koh_errors = c(
  89101431, # Empty koh string
  84228145, # Pretty sure this one is wrong
  47176946, # pretty sure this is wrong
  83427584,
  11662317,
  24397244, # strong case that this is wrong
  27832234, # another strong case,
  231530875,
  40660992,
  73401849,
  184999526,
  100132569,
  133647119
)

source("inst/make_exhuastive_test_cases.R")

repeat_vs_mh_ambiguity = c(109160237, 109711828)

deliberate_changes_from_koh = c(162229603) # This was bad design choice

devtools::load_all()

split_xx = test_indel_categorization()
s_diffs = which(split_xx$Koh_89 != split_xx$koh_orig_edited)
split_xx[s_diffs, ] -> diff_table
View(split_xx[s_diffs, ] |> dplyr::filter(!grepl(":M", Koh89.annotate.class)))

### Check the cases where we both call microhomology but different amounts
### subset of the following:
View(
  split_xx[s_diffs, ] |>
    dplyr::filter(ins_or_del == "d") |>
    dplyr::filter(grepl(":M", Koh_89))
)
# 0 cases

View(split_xx)

dplyr::filter(
  split_xx,
  ins_or_del == "d" &
    U_seq_count_in_indel_seq > 1 &
    R_outside_ins_or_del_seq >= U_seq_count_in_indel_seq
) -> foo


dplyr::filter(
  split_xx,
  ins_or_del == "d" &
    U_seq_count_in_indel_seq > 1 &
    R_outside_ins_or_del_seq < U_seq_count_in_indel_seq
) -> bar

dplyr::filter(
  split_xx,
  ins_or_del == "i" &
    U_seq_count_in_indel_seq > 1 &
    R_outside_ins_or_del_seq >= U_seq_count_in_indel_seq
) -> fii

View(fii)

dplyr::filter(
  split_xx,
  ins_or_del == "i" &
    U_seq_count_in_indel_seq > 1 &
    R_outside_ins_or_del_seq < U_seq_count_in_indel_seq
) -> bir

View(bir)

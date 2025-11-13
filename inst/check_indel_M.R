source("inst/make_exhuastive_test_cases.R")
library(dplyr)

devtools::load_all()
split_xx = test_indel_categorization()

filter(split_xx, ins_or_del == "d") -> split_xx

x476diffs = which(split_xx$Koh_476 != split_xx$Koh476.annotate.class)
length(x476diffs)
wdiff = split_xx[x476diffs, ] |>
  dplyr::select(
    -prev_COSMIC_83,
    -COSMIC_83,
    -koh_orig_edited,
    -Koh_89,
    -Koh89.annotate.class
  )
View(wdiff)

cosmh = wdiff
View(cosmh)

kohmh = wdiff
View(kohmh)


# View(dplyr::select(wdiff, Koh_476, Koh476.annotate.class))

###################################3

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

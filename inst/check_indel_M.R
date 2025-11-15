source("inst/make_exhuastive_test_cases.R")
library(dplyr)

tt = function(i_or_d) {
  devtools::load_all()
  split_xx = test_indel_categorization()

  filter(split_xx, ins_or_del == i_or_d) -> split_xx

  x476diffs = which(split_xx$Koh_476 != split_xx$Koh476.annotate.class)
  message("num diffs = ", length(x476diffs))
  wdiff = split_xx[x476diffs, ] |>
    dplyr::select(
      -prev_COSMIC_83,
      -COSMIC_83,
      -koh_orig_edited,
      -Koh_89,
      -Koh89.annotate.class
    )
  wdiff
}

deldif2 = tt("d")
insdif = tt("i")

err_in_koh = c(
  101303463, # micrhomology length is really 5, not 4
)

not_sure = c(100132569, 110914761)

frommo = read.csv("inst/annotations.unique.rows.csv", skip = 1)
sup7 = read.csv("inst/koh476types_from_sup_tab.csv")

##############################

# write.csv(koh476_from_sheet, "inst/koh476types_from_sup_tab.csv", row.names = F)
koh476_from_sheet = read.csv("inst/koh476types_from_sup_tab.csv", header = TRUE)
# koh476_from_sheet = koh476_from_sheet[, 2, drop = FALSE]

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

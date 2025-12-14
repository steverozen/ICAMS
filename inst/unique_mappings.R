vcf1 = read.csv("tests/testthat/testdata/long_test_vs_indelsiglib.csv")
vcf2 = read.csv(
  "tests/testthat/testdata/justified.mutations.ICAMS.to.debug.txt",
  sep = "\t"
)
cols = intersect(colnames(vcf1), colnames(vcf2))
vcf3 = rbind(vcf1[, cols], vcf2[, cols])

# vcf4 = read.csv("tests/testthat/testdata/M0_tests.tsv", sep = "\t")
# vcf5 = rbind(vcf3[, cols], vcf4[, cols])

retval1 = ICAMS:::categorize_indels_in_vcf(vcf3)
data.table::rbindlist(retval1, fill = TRUE) -> xx

length(unique(xx$Koh_476))

xx83 = xx |>
  dplyr::select(COSMIC_83, Koh_476) |>
  dplyr::distinct() |>
  group_by(Koh_476) |>
  mutate(count_83 = dplyr::n()) |>
  dplyr::filter(count_83 > 1)

xx89 = xx |>
  dplyr::select(Koh_89, Koh_476) |>
  dplyr::distinct() |>
  group_by(Koh_476) |>
  mutate(count_89 = dplyr::n()) |>
  dplyr::filter(count_89 > 1)

write.csv(xx89, "inst/multimap_to_89.csv")
write.csv(xx83, "inst/multimap_to_83.csv")
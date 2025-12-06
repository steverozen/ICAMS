test_that("test_long_test_vs_indelsiglib2", {
  # Warning: very slow
  vcf = read.csv(
    "testdata/justified.mutations.ICAMS.to.debug.txt",
    sep = "\t"
  ) |>
    dplyr::mutate(orig_icams_476 = Koh_476) |>
    dplyr::select(
      c(
        "Sample",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "seq.context.width",
        "seq.context",
        "pos_shift",
        "Koh_code_476",
        "orig_icams_476"
      )
    )
  rr = 1:nrow(vcf)
  retval1 = ICAMS:::categorize_indels_in_vcf(vcf[rr, ])

  cbind(vcf[rr, ], data.table::rbindlist(retval1, fill = TRUE)) -> xx
  qq = which(xx$Koh_476 != xx$Koh_code_476)
  if (length(qq) > 0) {
    View(dplyr::distinct(xx[qq, ] |> dplyr::select(Koh_476, Koh_code_476)))
  }
  expect_equal(length(qq), 0)
})

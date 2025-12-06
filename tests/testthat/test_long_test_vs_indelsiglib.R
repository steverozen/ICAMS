test_that("test_long_test_vs_indelsiglib", {
  mock_vcf = read.csv("testdata/long_test_vs_indelsiglib.csv")
  retval1 = ICAMS:::categorize_indels_in_vcf(mock_vcf)
  cbind(mock_vcf, data.table::rbindlist(retval1, fill = TRUE)) -> xx
  xx$koh_orig_edited = make_koh_open_intervals(xx)
  zz = which(xx$Koh_89 != xx$koh_orig_edited)
  qq = which(xx$Koh_476 != xx$Koh476.annotate.class)
  #     View(xx[zz, ] |> dplyr::select(Koh_89, koh_orig_edited))
  #     View(xx[qq, ] |> dplyr::select(Koh_476, Koh476.annotate.class))
  expect_equal(length(qq), 0)
  expect_equal(length(zz), 1)
  expect_snapshot(xx[zz, ])
})

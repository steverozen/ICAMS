test_that("test_long_test_vs_indelsiglib", {
  mock_vcf = read.csv("testdata/long_test_vs_indelsiglib.csv")
  retval1 = ICAMS:::categorize_indels_in_vcf(mock_vcf)
  cbind(mock_vcf, data.table::rbindlist(retval1, fill = TRUE)) -> xx
  xx$koh_orig_edited = make_koh_open_intervals(xx)
  zz = which(xx$Koh_89 != xx$koh_orig_edited)
  expect_snapshot(xx[zz, ])
})

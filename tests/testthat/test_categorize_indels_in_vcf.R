test_that("categorize_indels_in_vcf", {
  fixture_path <- testthat::test_path(
    "testdata",
    "categorize_many_indels_test_input.csv"
  )
  dd = read.csv(
    fixture_path # "testdata/categorize_many_indels_test_input.csv"
  )

  yy = data.table::rbindlist(ICAMS:::categorize_indels_in_vcf(dd), fill = TRUE)

  expect_snapshot(t(yy))
})

test_that("categorize_indels_in_vcf", {
  dd = read.csv(
    "testdata/categorize_many_indels_test_input.csv"
  )

  yy = data.table::rbindlist(ICAMS:::categorize_indels_in_vcf(dd), fill = TRUE)

  expect_snapshot(yy)
})

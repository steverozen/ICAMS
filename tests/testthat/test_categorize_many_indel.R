test_that("categorize_many_indels", {
  source("../../inst/xcategorize_1_justified_indel.R")
  dd = read.csv(
    "testdata/categorize_many_indels_test_input.csv"
  )

  yy = data.table::rbindlist(ICAMS:::categorize_many_indels(dd), fill = TRUE)

  expect_snapshot(yy)
})

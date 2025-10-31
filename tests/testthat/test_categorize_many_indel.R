test_that("categorize_many_indels", {
  dd = read.csv(
    "testdata/categorize_many_indels_test_input.csv"
  )

  yy = data.table::rbindlist(ICAMS:::categorize_many_indels(dd), fill = TRUE)

  expect_snapshot(yy)
})

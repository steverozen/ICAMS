context("PlotCat192")

test_that("PlotCat192 function is working properly", {
  par(mar = c(6, 4, 4, 1))
  catalog <- ReadCat192("testdata/regress.cat.192.csv")
  out <- PlotCat192(catalog[, 1, drop = FALSE], "test")
  expect_equal(out, TRUE)
})

context("PlotCat192Strand")

test_that("PlotCat192Strand function is working properly", {
  par(mar = c(5, 8, 5, 1))
  catalog <- ReadCat192("testdata/regress.cat.192.csv")
  out <- PlotCat192Strand(catalog[, 1, drop = FALSE], "test")
  expect_equal(out, TRUE)
})

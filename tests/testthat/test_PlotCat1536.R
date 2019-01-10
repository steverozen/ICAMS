context("PlotCat1536")

test_that("PlotCat1536 function is working properly", {
  catalog <- ReadCat1536("testdata/regress.cat.1536.csv")
  out <-
    PlotCat1536(catalog[, 1, drop = FALSE], "test", abundance = .abundance.5bp)
  expect_equal(out, TRUE)
})

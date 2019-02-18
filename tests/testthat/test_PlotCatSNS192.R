context("PlotCatSNS192")

test_that("PlotCatSNS192 function is working properly", {
  par(mar = c(6, 4, 4, 1))
  catalog <- ReadCatSNS192("testdata/regress.cat.sns.192.csv")
  out <- PlotCatSNS192(catalog[, 1, drop = FALSE], "test")
  expect_equal(out, TRUE)
})

context("PlotCatSNS192Strand")

test_that("PlotCatSNS192Strand function is working properly", {
  par(mar = c(5, 8, 5, 1))
  catalog <- ReadCatSNS192("testdata/regress.cat.sns.192.csv")
  out <- PlotCatSNS192Strand(catalog[, 1, drop = FALSE], "test")
  expect_equal(out, TRUE)
})

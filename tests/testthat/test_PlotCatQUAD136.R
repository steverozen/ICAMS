context("PlotCatQUAD136")

test_that("PlotCatQUAD136 function is working properly", {
  catalog <- ReadCatQUAD136("testdata/regress.cat.quad.136.csv")
  out <- PlotCatQUAD136(catalog[, 1, drop = FALSE], id = "test",
                      type = "density", abundance = .abundance.4bp)
  expect_equal(out, TRUE)
})

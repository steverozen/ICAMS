context("PlotCat96")

test_that("PlotCat96 function is working properly", {
  par(mar = c(11.5, 5.1, 4.6, 3))
  catalog <- ReadCat96("testdata/regress.cat.96.csv")
  out <- PlotCat96(catalog[, 1, drop = FALSE], id = "test",
            type = "density", abundance = .abundance.3bp)
  expect_equal(out, TRUE)
})

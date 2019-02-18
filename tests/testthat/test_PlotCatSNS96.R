context("PlotCatSNS96")

test_that("PlotCatSNS96 function is working properly", {
  par(mar = c(11.5, 5.1, 4.6, 3))
  catalog <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")
  out <- PlotCatSNS96(catalog[, 1, drop = FALSE], id = "test",
                   type = "density", grid = TRUE, abundance = .abundance.3bp.GRCh37)
  expect_equal(out, TRUE)
})

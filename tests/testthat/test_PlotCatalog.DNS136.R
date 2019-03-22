context("PlotCatalog.DNS136")

test_that("PlotCatalog.DNS136 function is working properly", {
  par(oma = c(2, 2, 2, 0))
  catalog <- ReadCatalog("testdata/regress.cat.dns.136.csv",
                         ref.genome = "GRCh37",
                         region = "genome", type = "counts")
  catalog$catalog <- catalog$catalog[, 1, drop = FALSE]
  out <- PlotCatalog(catalog)
  expect_equal(out, TRUE)
})

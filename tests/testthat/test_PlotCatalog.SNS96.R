context("PlotCatalog.SNS96")

test_that("PlotCatalog.SNS96 function is working properly", {
  par(mar = c(11.5, 5.1, 4.6, 3))
  catalog <- ReadCatalog("testdata/regress.cat.sns.96.csv", ref.genome = "GRCh37",
                         region = "genome", type = "counts")
  catalog$catalog <- catalog$catalog[, 1, drop = FALSE]
  out <- PlotCatalog(catalog)
  expect_equal(out, TRUE)
})

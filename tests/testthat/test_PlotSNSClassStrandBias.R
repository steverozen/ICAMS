context("PlotSNSClassStrandBias")

test_that("PlotSNSClassStrandBias function is working properly", {
  par(mar = c(5, 8, 5, 1))
  catalog <- ReadCatalog("testdata/regress.cat.sns.192.csv", ref.genome = "GRCh37",
                           region = "genome", type = "counts")
  catalog$catalog <- catalog$catalog[, 1, drop = FALSE]
  out <- PlotSNSClassStrandBias(catalog)
  expect_equal(out, TRUE)
})

context("PlotCatalog.DNS78")

test_that("PlotCatalog.DNS78 function is working properly", {
  par(mar = c(6, 4, 6, 2))
  catalog <- ReadCatalog("testdata/regress.cat.dns.78.csv",
                          ref.genome = "GRCh37",
                          region = "genome", type = "counts")
  cat <- catalog[, 1, drop = FALSE]
  cat <- PreserveCatalogAttribute(catalog, cat)
  out <- PlotCatalog(cat)
  expect_equal(out, TRUE)
})

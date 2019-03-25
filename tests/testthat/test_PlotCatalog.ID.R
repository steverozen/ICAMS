context("PlotCatalog.ID")

test_that("PlotCatalog.ID function is working properly", {
  par(mar = c(7, 4, 7, 3))
  catalog <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                       ref.genome = "GRCh37",
                       region = "genome", type = "counts")
  cat <- catalog[, 1, drop = FALSE]
  cat <- PreserveCatalogAttribute(catalog, cat)
  out <- PlotCatalog(cat)
  expect_equal(out, TRUE)
})

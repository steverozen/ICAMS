context("PlotCatalog.SNS96")

test_that("PlotCatalog.SNS96 function is working properly", {
  par(mar = c(8, 5, 5, 1))
  catalog <- ReadCatalog("testdata/regress.cat.sns.96.csv",
                         ref.genome = "GRCh37",
                         region = "genome", catalog.type = "counts")
  cat <- catalog[, 1, drop = FALSE]
  cat <- PreserveCatalogAttribute(catalog, cat)
  out <- PlotCatalog(cat)
  expect_equal(out, TRUE)
})

context("PlotCatalog.SNS192")

test_that("PlotCatalog.SNS192 function is working properly", {
  par(mar = c(6, 4, 4, 1))
  catalog <- ReadCatalog("testdata/regress.cat.sns.192.csv", ref.genome = "GRCh37",
                         region = "genome", catalog.type = "counts")
  cat <- catalog[, 1, drop = FALSE]
  cat <- PreserveCatalogAttribute(catalog, cat)
  out <- PlotCatalog(cat)
  expect_equal(out, TRUE)
})

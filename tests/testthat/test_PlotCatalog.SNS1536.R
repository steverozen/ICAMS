context("PlotCatalog.SNS1536")

test_that("PlotCatalog.SNS1536 function is working properly", {
  catalog <- ReadCatalog("testdata/regress.cat.sns.1536.csv", ref.genome = "GRCh37",
                            region = "genome", type = "counts")
  cat <- catalog[, 1, drop = FALSE]
  cat <- PreserveCatalogAttribute(catalog, cat)
  out <- PlotCatalog(cat)
  expect_equal(out, TRUE)
})

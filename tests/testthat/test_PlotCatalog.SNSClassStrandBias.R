context("PlotCatalog.SNSClassStrandBias")

test_that("PlotCatalog.SNSClassStrandBias function is working properly", {
  par(mar = c(5, 8, 5, 1))
  catalog <- ReadCatalog("testdata/regress.cat.sns.192.csv", ref.genome = "GRCh37",
                           region = "genome", catalog.type = "counts")
  cat <- catalog[, 1, drop = FALSE]
  cat <- PreserveCatalogAttribute(catalog, cat)
  out <- PlotCatalog(cat, strandbias = TRUE)
  expect_equal(out, TRUE)
})

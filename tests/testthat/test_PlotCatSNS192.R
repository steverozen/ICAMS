context("PlotCatSNS192")

test_that("PlotCatSNS192 function is working properly", {
  par(mar = c(6, 4, 4, 1))
  catalog <- ReadCatSNS192("testdata/regress.cat.sns.192.csv")
  out <- PlotCatSNS192(catalog[, 1, drop = FALSE], type = "counts", id = "test")
  expect_equal(out, TRUE)

  signature.catalog <-
    TransformCatalog(catalog, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome",
                     which.n = 3, source.type = "counts",
                     target.type = "signature")
  out1 <- PlotCatSNS192(signature.catalog[, 1, drop = FALSE],
                        type = "signature", id = "test")
  expect_equal(out1, TRUE)

})

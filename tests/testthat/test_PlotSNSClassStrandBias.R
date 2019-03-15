context("PlotSNSClassStrandBias")

test_that("PlotSNSClassStrandBias function is working properly", {
  par(mar = c(5, 8, 5, 1))
  catalog <- ReadCatSNS192("testdata/regress.cat.sns.192.csv")
  out <- PlotSNSClassStrandBias(catalog[, 1, drop = FALSE],
                                type = "counts", id = "test")
  expect_equal(out, TRUE)

  signature.catalog <-
    TransformCatalog(catalog, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome",
                     which.n = 3, source.type = "counts",
                     target.type = "signature")
  out1 <-
    PlotSNSClassStrandBias(signature.catalog[, 1, drop = FALSE],
                           id = "test", type = "signature")
  expect_equal(out1, TRUE)
})

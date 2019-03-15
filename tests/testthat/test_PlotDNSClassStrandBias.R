context("PlotDNSClassStrandBias")

test_that("PlotDNSClassStrandBias function is working properly", {
  par(mar = c(5, 8, 5, 1))
  catalog <- ReadCatDNS144("testdata/regress.cat.dns.144.csv")
  out <- PlotDNSClassStrandBias(catalog[, 1, drop = FALSE],
                                type = "counts", id = "test")
  expect_equal(out, TRUE)

  signature.catalog <-
    TransformCatalog(catalog, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome",
                     which.n = 2, source.type = "counts",
                     target.type = "signature")
  out1 <-
    PlotDNSClassStrandBias(signature.catalog[, 1, drop = FALSE],
                           id = "test", type = "signature")
  expect_equal(out1, TRUE)
})

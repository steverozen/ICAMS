context("PlotCatDNS78")

test_that("PlotCatDNS78 function is working properly", {
  par(mar = c(6, 4, 6, 2))
  catalog <- ReadCatDNS78("testdata/regress.cat.dns.78.csv")
  out <-
    PlotCatDNS78(catalog[, 1, drop = FALSE],
                 id = "test",
                 type = "counts")
  expect_equal(out, TRUE)

  density.catalog <-
    TransformCatalog(catalog, source.abundance = "GRCh37.genome",
                     which.n = 2, source.type = "counts", target.type = "density")
  out1 <-
    PlotCatDNS78(density.catalog[, 1, drop = FALSE],
                 id = "test", type = "density")
  expect_equal(out1, TRUE)

  signature.catalog <-
    TransformCatalog(catalog, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome",
                     which.n = 2, source.type = "counts", target.type = "signature")
  out2 <-
    PlotCatDNS78(signature.catalog[, 1, drop = FALSE],
                 id = "test", type = "signature")
  expect_equal(out2, TRUE)
})

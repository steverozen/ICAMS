context("PlotCatDNS136")

test_that("PlotCatDNS136 function is working properly", {
  par(oma = c(2, 2, 2, 0))
  catalog <- ReadCatDNS136("testdata/regress.cat.dns.136.csv")
  out <-
    PlotCatDNS136(catalog[, 1, drop = FALSE],
                  id = "test",
                  type = "counts")
  expect_equal(out, TRUE)

  density.catalog <-
    TransformCatalog(catalog, source.abundance = "GRCh37.genome",
                     which.n = 4, source.type = "counts",
                     target.type = "density")
  out1 <-
    PlotCatDNS136(density.catalog[, 1, drop = FALSE],
                 id = "test", type = "density")
  expect_equal(out1, TRUE)
})

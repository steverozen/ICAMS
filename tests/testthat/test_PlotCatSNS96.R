context("PlotCatSNS96")

test_that("PlotCatSNS96 function is working properly", {
  par(mar = c(11.5, 5.1, 4.6, 3))
  catalog <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")
  out <-
    PlotCatSNS96(catalog[, 1, drop = FALSE],
                 id = "test", type = "counts")
  expect_equal(out, TRUE)

  density.catalog <-
    TransformCatalog(catalog, source.abundance = "GRCh37.genome",
                     which.n = 3, source.type = "counts", target.type = "density")
  out1 <-
    PlotCatSNS96(density.catalog[, 1, drop = FALSE],
                 id = "test", type = "density")
  expect_equal(out1, TRUE)

  signature.catalog <-
    TransformCatalog(catalog, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome",
                     which.n = 3, source.type = "counts", target.type = "signature")
  out2 <-
    PlotCatSNS96(signature.catalog[, 1, drop = FALSE],
                 id = "test", type = "signature")
  expect_equal(out2, TRUE)
})

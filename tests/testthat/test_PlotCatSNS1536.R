context("PlotCatSNS1536")

test_that("PlotCatSNS1536 function is working properly", {
  catalog <- ReadCatSNS1536("testdata/regress.cat.sns.1536.csv")
  out <-
    PlotCatSNS1536(catalog[, 1, drop = FALSE], type = "counts",
                   id = "test")
  expect_equal(out, TRUE)

  density.catalog <-
    TransformCatalog(catalog, source.abundance = "GRCh37.genome",
                     which.n = 5, source.type = "counts", target.type = "density")
  out1 <-
    PlotCatSNS1536(density.catalog[, 1, drop = FALSE],
                 id = "test", type = "density")
  expect_equal(out1, TRUE)

  signature.catalog <-
    TransformCatalog(catalog, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome",
                     which.n = 5, source.type = "counts", target.type = "signature")
  out2 <-
    PlotCatSNS1536(signature.catalog[, 1, drop = FALSE],
                 id = "test", type = "signature")
  expect_equal(out2, TRUE)
})

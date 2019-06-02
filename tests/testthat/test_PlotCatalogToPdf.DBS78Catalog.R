context("PlotCatalogToPdf.DBS78Catalog")

test_that("PlotCatalogToPdf.DBS78Catalog function is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.78.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog.counts, file = "PlotCatDBS78.counts.test.pdf")
  expect_equal(out, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density, file = "PlotCatDBS78.density.test.pdf")
  expect_equal(out, TRUE)

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = "PlotCatDBS78.counts.signature.test.pdf")
  expect_equal(out, TRUE)

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     file = "PlotCatDBS78.density.signature.test.pdf")
  expect_equal(out, TRUE)

  unlink("PlotCatDBS78.counts.test.pdf")
  unlink("PlotCatDBS78.density.test.pdf")
  unlink("PlotCatDBS78.counts.signature.test.pdf")
  unlink("PlotCatDBS78.density.signature.test.pdf")
  graphics.off()
})
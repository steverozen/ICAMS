context("PlotCatalogToPdf.DBS144Catalog")

test_that("PlotCatalogToPdf.DBS144Catalog function is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.144.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog.counts,
                     file = "PlotDBS144Catalog.counts.test.pdf")
  expect_equal(out, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density,
                     file = "PlotDBS144Catalog.density.test.pdf")
  expect_equal(out, TRUE)

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = "PlotDBS144Catalog.counts.signature.test.pdf")
  expect_equal(out, TRUE)

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     file = "PlotDBS144Catalog.density.signature.test.pdf")
  expect_equal(out, TRUE)

  unlink("PlotDBS144Catalog.counts.test.pdf")
  unlink("PlotDBS144Catalog.density.test.pdf")
  unlink("PlotDBS144Catalog.counts.signature.test.pdf")
  unlink("PlotDBS144Catalog.density.signature.test.pdf")
  graphics.off()
})
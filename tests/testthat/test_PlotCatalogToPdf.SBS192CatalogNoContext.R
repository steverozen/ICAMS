context("PlotCatalogToPdf.SBS192CatalogNoContext")

test_that("PlotCatalogToPdf.SBS192CatalogNoContext function is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.192.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog.counts,
                     file = "PlotSBS192CatalogNoContext.counts.test.pdf",
                     no.context = TRUE)
  expect_equal(out, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density,
                     file = "PlotSBS192CatalogNoContext.density.test.pdf",
                     no.context = TRUE)
  expect_equal(out, TRUE)

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = "PlotSBS192CatalogNoContext.counts.signature.test.pdf",
                     no.context = TRUE)
  expect_equal(out, TRUE)

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     file = "PlotSBS192CatalogNoContext.density.signature.test.pdf",
                     no.context = TRUE)
  expect_equal(out, TRUE)

  unlink("PlotSBS192CatalogNoContext.counts.test.pdf")
  unlink("PlotSBS192CatalogNoContext.density.test.pdf")
  unlink("PlotSBS192CatalogNoContext.counts.signature.test.pdf")
  unlink("PlotSBS192CatalogNoContext.density.signature.test.pdf")
  graphics.off()
})

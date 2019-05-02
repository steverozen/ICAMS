context("PlotCatalogToPdf.DNS78")

test_that("PlotCatalogToPdf.DNS78 function is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dns.78.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog.counts, filename = "PlotCatDNS78.counts.test.pdf")
  expect_equal(out, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density, filename = "PlotCatDNS78.density.test.pdf")
  expect_equal(out, TRUE)

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     filename = "PlotCatDNS78.counts.signature.test.pdf")
  expect_equal(out, TRUE)

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     filename = "PlotCatDNS78.density.signature.test.pdf")
  expect_equal(out, TRUE)

  unlink("PlotCatDNS78.counts.test.pdf")
  unlink("PlotCatDNS78.density.test.pdf")
  unlink("PlotCatDNS78.counts.signature.test.pdf")
  unlink("PlotCatDNS78.density.signature.test.pdf")
  graphics.off()
})

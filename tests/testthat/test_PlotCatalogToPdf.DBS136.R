context("PlotCatalogToPdf.DNS136")

test_that("PlotCatalogToPdf.DNS136 function is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dns.136.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog.counts, file = "PlotCatDNS136.counts.test.pdf")
  expect_equal(out, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density, file = "PlotCatDNS136.density.test.pdf")
  expect_equal(out, TRUE)

  unlink("PlotCatDNS136.counts.test.pdf")
  unlink("PlotCatDNS136.density.test.pdf")
  graphics.off()
})

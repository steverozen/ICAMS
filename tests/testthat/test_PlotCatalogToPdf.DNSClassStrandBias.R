context("PlotCatalogToPdf.DNSClassStrandBias")

test_that("PlotCatalogToPdf.DNSClassStrandBias function is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dns.144.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog.counts, filename = "PlotDNSClassStrandBias.counts.test.pdf",
                     strandbias = TRUE)
  expect_equal(out, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density, filename = "PlotDNSClassStrandBias.density.test.pdf",
                     strandbias = TRUE)
  expect_equal(out, TRUE)

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     filename = "PlotDNSClassStrandBias.counts.signature.test.pdf",
                     strandbias = TRUE)
  expect_equal(out, TRUE)

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     filename = "PlotDNSClassStrandBias.density.signature.test.pdf",
                     strandbias = TRUE)
  expect_equal(out, TRUE)

  unlink("PlotDNSClassStrandBias.counts.test.pdf")
  unlink("PlotDNSClassStrandBias.density.test.pdf")
  unlink("PlotDNSClassStrandBias.counts.signature.test.pdf")
  unlink("PlotDNSClassStrandBias.density.signature.test.pdf")
  graphics.off()
})

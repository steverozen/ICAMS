context("PlotCatalogToPdf.SNSClassStrandBias")

test_that("PlotCatalogToPdf.SNSClassStrandBias function is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sns.192.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog.counts,
                     filename = "PlotSNSClassStrandBias.counts.test.pdf",
                     strandbias = TRUE)
  expect_equal(out, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density,
                     filename = "PlotSNSClassStrandBias.density.test.pdf",
                     strandbias = TRUE)
  expect_equal(out, TRUE)

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     filename = "PlotSNSClassStrandBias.counts.signature.test.pdf",
                     strandbias = TRUE)
  expect_equal(out, TRUE)

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     filename = "PlotSNSClassStrandBias.density.signature.test.pdf",
                     strandbias = TRUE)
  expect_equal(out, TRUE)

  unlink("PlotSNSClassStrandBias.counts.test.pdf")
  unlink("PlotSNSClassStrandBias.density.test.pdf")
  unlink("PlotSNSClassStrandBias.counts.signature.test.pdf")
  unlink("PlotSNSClassStrandBias.density.signature.test.pdf")
  graphics.off()
})

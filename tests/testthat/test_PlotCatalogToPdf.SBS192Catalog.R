context("PlotCatalogToPdf.SBS192Catalog")

test_that("PlotCatalogToPdf.SBS192Catalog function is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.192.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog.counts, 
                     file = paste0(tempdir(), "\\PlotCatSBS192.counts.test.pdf"))
  out1 <-
    PlotCatalogToPdf(catalog.counts,
                     file = paste0(tempdir(), "\\PlotCatSBS192NoContext.counts.test.pdf"),
                     no.context = TRUE)
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density, file = paste0(tempdir(), "\\PlotCatSBS192.density.test.pdf"))
  out1 <-
    PlotCatalogToPdf(catalog.density,
                     file = paste0(tempdir(), "\\PlotCatSBS192NoContext.density.test.pdf"),
                     no.context = TRUE)
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = paste0(tempdir(), "\\PlotCatSBS192.counts.signature.test.pdf"))
  out1 <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = paste0(tempdir(), "\\PlotCatSBS192NoContext.counts.signature.test.pdf"),
                     no.context = TRUE)
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     file = paste0(tempdir(), "\\PlotCatSBS192.density.signature.test.pdf"))
  out1 <-
    PlotCatalogToPdf(catalog.density.signature,
                     file = paste0(tempdir(), "\\PlotCatSBS192NoContext.density.signature.test.pdf"),
                     no.context = TRUE)
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)

  unlink(paste0(tempdir(), "\\PlotCatSBS192.counts.test.pdf"))
  unlink(paste0(tempdir(), "\\PlotCatSBS192.density.test.pdf"))
  unlink(paste0(tempdir(), "\\PlotCatSBS192.counts.signature.test.pdf"))
  unlink(paste0(tempdir(), "\\PlotCatSBS192.density.signature.test.pdf"))
  unlink(paste0(tempdir(), "\\PlotCatSBS192NoContext.counts.test.pdf"))
  unlink(paste0(tempdir(), "\\PlotCatSBS192NoContext.density.test.pdf"))
  unlink(paste0(tempdir(), "\\PlotCatSBS192NoContext.counts.signature.test.pdf"))
  unlink(paste0(tempdir(), "\\PlotCatSBS192NoContext.density.signature.test.pdf"))
  graphics.off()
})

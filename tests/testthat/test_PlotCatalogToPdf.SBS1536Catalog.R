context("PlotCatalogToPdf.SBS1536Catalog")

test_that("PlotCatalogToPdf.SBS1536Catalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.1536.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog.counts, 
                     file = file.path(tempdir(), "PlotCatSBS1536.counts.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density, file = file.path(tempdir(), "PlotCatSBS1536.density.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = file.path(tempdir(), "PlotCatSBS1536.counts.signature.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     file = file.path(tempdir(), "PlotCatSBS1536.density.signature.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  unlink(file.path(tempdir(), "PlotCatSBS1536.counts.test.pdf"))
  unlink(file.path(tempdir(), "PlotCatSBS1536.density.test.pdf"))
  unlink(file.path(tempdir(), "PlotCatSBS1536.counts.signature.test.pdf"))
  unlink(file.path(tempdir(), "PlotCatSBS1536.density.signature.test.pdf"))
  graphics.off()
})

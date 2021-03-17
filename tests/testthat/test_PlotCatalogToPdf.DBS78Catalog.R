context("PlotCatalogToPdf.DBS78Catalog")

test_that("PlotCatalogToPdf.DBS78Catalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/BTSG_WGS_PCAWG.dbs.78.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  out <-
    PlotCatalogToPdf(catalog.counts, 
                     file = file.path(tempdir(), "PlotCatDBS78.counts.test.pdf"))
  out1 <-
    PlotCatalogToPdf(catalog.counts, 
                     file = file.path(tempdir(), "PlotCatDBS78.counts.test2.pdf"),
                     grid = FALSE, upper = FALSE, xlabels = FALSE, ylabels = FALSE)
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density, 
                     file = file.path(tempdir(), "PlotCatDBS78.density.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = file.path(tempdir(), "PlotCatDBS78.counts.signature.test.pdf"))
  colnames(catalog.counts.signature) <- rep("", ncol(catalog.counts.signature))
  out1 <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = file.path(tempdir(), "PlotCatDBS78.counts.signature.test2.pdf"),
                     grid = FALSE, upper = FALSE, xlabels = FALSE, ylabels = FALSE)
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     file = file.path(tempdir(), "PlotCatDBS78.density.signature.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  unlink(file.path(tempdir(), "PlotCatDBS78.counts.test.pdf"))
  unlink(file.path(tempdir(), "PlotCatDBS78.counts.test2.pdf"))
  unlink(file.path(tempdir(), "PlotCatDBS78.density.test.pdf"))
  unlink(file.path(tempdir(), "PlotCatDBS78.counts.signature.test.pdf"))
  unlink(file.path(tempdir(), "PlotCatDBS78.counts.signature.test2.pdf"))
  unlink(file.path(tempdir(), "PlotCatDBS78.density.signature.test.pdf"))
  graphics.off()
})

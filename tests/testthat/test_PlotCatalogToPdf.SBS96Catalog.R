context("PlotCatalogToPdf.SBS96Catalog")

test_that("PlotCatalogToPdf.SBS96Catalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/BTSG_WGS_PCAWG.sbs.96.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  out <-
    PlotCatalogToPdf(catalog.counts, 
                     file = file.path(tempdir(), "PlotCatSBS96.counts.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density, file = file.path(tempdir(), "PlotCatSBS96.density.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = file.path(tempdir(), "PlotCatSBS96.counts.signature.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     file = file.path(tempdir(), "PlotCatSBS96.density.signature.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  if (Sys.getenv("ICAMS.SAVE.TEST.PDF") != "") {
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatSBS96.counts.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS96",
                                 "PlotCatSBS96.counts.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatSBS96.density.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS96",
                                 "PlotCatSBS96.density.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatSBS96.counts.signature.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS96",
                                 "PlotCatSBS96.counts.signature.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatSBS96.density.signature.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS96",
                                 "PlotCatSBS96.density.signature.test.pdf"))
  } else {
    unlink(file.path(tempdir(), "PlotCatSBS96.counts.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatSBS96.density.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatSBS96.counts.signature.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatSBS96.density.signature.test.pdf"))
  }
  graphics.off()
})

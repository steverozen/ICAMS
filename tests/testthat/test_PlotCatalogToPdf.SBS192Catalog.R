context("PlotCatalogToPdf.SBS192Catalog")

test_that("PlotCatalogToPdf.SBS192Catalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/BTSG_WGS_PCAWG.sbs.192.csv",
                                ref.genome = "GRCh37",
                                region = "transcript", catalog.type = "counts")
  out <-
    PlotCatalogToPdf(catalog.counts, 
                     file = file.path(tempdir(), "PlotCatSBS192.counts.test.pdf"))
  out1 <-
    PlotCatalogToPdf(catalog.counts,
                     file = file.path(tempdir(), "PlotCatalog.SBS12.counts.test.pdf"),
                     plot.SBS12 = TRUE)
  expect_equal(out$plot.success, TRUE) 
  expect_equal(out1$plot.success, TRUE) 

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density, file = file.path(tempdir(), "PlotCatSBS192.density.test.pdf"))
  out1 <-
    PlotCatalogToPdf(catalog.density,
                     file = file.path(tempdir(), "PlotCatalog.SBS12.density.test.pdf"),
                     plot.SBS12 = TRUE)
  expect_equal(out$plot.success, TRUE) 
  expect_equal(out1$plot.success, TRUE) 

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = file.path(tempdir(), "PlotCatSBS192.counts.signature.test.pdf"))
  out1 <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = file.path(tempdir(), "PlotCatalog.SBS12.counts.signature.test.pdf"),
                     plot.SBS12 = TRUE)
  expect_equal(out$plot.success, TRUE) 
  expect_equal(out1$plot.success, TRUE) 

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     file = file.path(tempdir(), "PlotCatSBS192.density.signature.test.pdf"))
  out1 <-
    PlotCatalogToPdf(catalog.density.signature,
                     file = file.path(tempdir(), "PlotCatalog.SBS12.density.signature.test.pdf"),
                     plot.SBS12 = TRUE)
  expect_equal(out$plot.success, TRUE) 
  expect_equal(out1$plot.success, TRUE) 
  
  if (Sys.getenv("ICAMS.SAVE.TEST.PDF") != "") {
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatSBS192.counts.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS192",
                                 "PlotCatSBS192.counts.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatSBS192.density.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS192",
                                 "PlotCatSBS192.density.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatSBS192.counts.signature.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS192",
                                 "PlotCatSBS192.counts.signature.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatSBS192.density.signature.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS192",
                                 "PlotCatSBS192.density.signature.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatalog.SBS12.counts.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS192",
                                 "PlotCatalog.SBS12.counts.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatalog.SBS12.density.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS192",
                                 "PlotCatalog.SBS12.density.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatalog.SBS12.counts.signature.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS192",
                                 "PlotCatalog.SBS12.counts.signature.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatalog.SBS12.density.signature.test.pdf"),
                to   = file.path("pdfs-for-comparision-SBS192",
                                 "PlotCatalog.SBS12.density.signature.test.pdf"))
  } else {
    unlink(file.path(tempdir(), "PlotCatSBS192.counts.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatSBS192.density.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatSBS192.counts.signature.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatSBS192.density.signature.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatalog.SBS12.counts.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatalog.SBS12.density.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatalog.SBS12.counts.signature.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatalog.SBS12.density.signature.test.pdf"))
  }
  graphics.off()
})

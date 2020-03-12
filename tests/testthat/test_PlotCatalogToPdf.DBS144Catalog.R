context("PlotCatalogToPdf.DBS144Catalog")

test_that("PlotCatalogToPdf.DBS144Catalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.144.csv",
                                ref.genome = "GRCh37",
                                region = "transcript",
                                catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog.counts,
                     file = file.path(tempdir(), "PlotCatDBS144.counts.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density,
                     file = file.path(tempdir(), "PlotCatDBS144.density.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  catalog.counts.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = file.path(tempdir(), "PlotCatDBS144.counts.signature.test.pdf"))
  expect_equal(out$plot.success, TRUE)

  catalog.density.signature <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "density.signature")
  out <-
    PlotCatalogToPdf(catalog.density.signature,
                     file = file.path(tempdir(), "PlotCatDBS144.density.signature.test.pdf"))
  expect_equal(out$plot.success, TRUE)
  
  if (Sys.getenv("ICAMS.SAVE.TEST.PDF") != "") {
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatDBS144.counts.test.pdf"),
                to   = file.path("pdfs-for-comparision-DBS144",
                                 "PlotCatDBS144.counts.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatDBS144.density.test.pdf"),
                to   = file.path("pdfs-for-comparision-DBS144",
                                 "PlotCatDBS144.density.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatDBS144.counts.signature.test.pdf"),
                to   = file.path("pdfs-for-comparision-DBS144",
                                 "PlotCatDBS144.counts.signature.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatDBS144.density.signature.test.pdf"),
                to   = file.path("pdfs-for-comparision-DBS144",
                                 "PlotCatDBS144.density.signature.test.pdf"))
  } else {
    unlink(file.path(tempdir(), "PlotCatDBS144.counts.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatDBS144.density.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatDBS144.counts.signature.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatDBS144.density.signature.test.pdf"))
  }
  
  graphics.off()
})

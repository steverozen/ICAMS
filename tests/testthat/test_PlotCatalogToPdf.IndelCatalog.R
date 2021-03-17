context("PlotCatalogToPdf.IndelCatalog")

test_that("PlotCatalogToPdf.IndelCatalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                       ref.genome = "GRCh37",
                       region = "genome", catalog.type = "counts")
  out <- PlotCatalogToPdf(catalog.counts, 
                          file = file.path(tempdir(), "PlotCatID.counts.test.pdf"))
  out1 <- PlotCatalogToPdf(catalog.counts, 
                           file = file.path(tempdir(), "PlotCatID.counts.test2.pdf"),
                           upper = FALSE, grid = FALSE, xlabels = FALSE, ylabels = FALSE)
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)

  catalog.counts.signature <-
    apply(catalog.counts, MARGIN = 2, function(x) x / sum(x))
  catalog.counts.signature <-
    as.catalog(catalog.counts.signature, ref.genome = "GRCh37",
               region = "genome", catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = file.path(tempdir(), "PlotCatID.counts.signature.test.pdf"))
  colnames(catalog.counts.signature) <- rep("", ncol(catalog.counts.signature))
  out1 <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = file.path(tempdir(), "PlotCatID.counts.signature.test2.pdf"),
                     grid = FALSE, upper = FALSE, xlabels = FALSE, ylabels = FALSE)
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)
  
  if (Sys.getenv("ICAMS.SAVE.TEST.PDF") != "") {
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatID.counts.test.pdf"),
                to   = file.path("pdfs-for-comparision-IndelCatalog",
                                 "PlotCatID.counts.test.pdf"))
    file.rename(from = file.path(tempdir(), 
                                 "PlotCatID.counts.signature.test.pdf"),
                to   = file.path("pdfs-for-comparision-IndelCatalog",
                                 "PlotCatID.counts.signature.test.pdf"))
  } else {
    unlink(file.path(tempdir(), "PlotCatID.counts.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatID.counts.test2.pdf"))
    unlink(file.path(tempdir(), "PlotCatID.counts.signature.test.pdf"))
    unlink(file.path(tempdir(), "PlotCatID.counts.signature.test2.pdf"))
  }
  graphics.off()
})

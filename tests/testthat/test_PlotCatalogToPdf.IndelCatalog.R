context("PlotCatalogToPdf.IndelCatalog")

test_that("PlotCatalogToPdf.IndelCatalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                       ref.genome = "GRCh37",
                       region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("Biliary-AdenoCA", 1 : 35)
  out <- PlotCatalogToPdf(catalog.counts, 
                          file = file.path(tempdir(), "PlotCatID.counts.test.pdf"))
  expect_equal(out, TRUE)

  catalog.counts.signature <-
    apply(catalog.counts, MARGIN = 2, function(x) x / sum(x))
  catalog.counts.signature <-
    as.catalog(catalog.counts.signature, ref.genome = "GRCh37",
               region = "genome", catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = file.path(tempdir(), "PlotCatID.counts.signature.test.pdf"))
  expect_equal(out, TRUE)

  unlink(file.path(tempdir(), "PlotCatID.counts.test.pdf"))
  unlink(file.path(tempdir(), "PlotCatID.counts.signature.test.pdf"))
  graphics.off()
})

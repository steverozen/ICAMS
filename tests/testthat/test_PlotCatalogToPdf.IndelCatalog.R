context("PlotCatalogToPdf.IndelCatalog")

test_that("PlotCatalogToPdf.IndelCatalog function is working properly", {
  catalog.counts <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                       ref.genome = "GRCh37",
                       region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("Biliary-AdenoCA", 1 : 35)
  out <- PlotCatalogToPdf(catalog.counts, file = "PlotCatID.test.pdf")
  expect_equal(out, TRUE)

  catalog.counts.signature <-
    apply(catalog.counts, MARGIN = 2, function(x) x / sum(x))
  catalog.counts.signature <-
    as.catalog(catalog.counts.signature, ref.genome = "GRCh37",
               region = "genome", catalog.type = "counts.signature")
  out <-
    PlotCatalogToPdf(catalog.counts.signature,
                     file = "PlotCatID.counts.signature.test.pdf")
  expect_equal(out, TRUE)

  unlink("PlotCatID.test.pdf")
  unlink("PlotCatID.counts.signature.test.pdf")
  graphics.off()
})

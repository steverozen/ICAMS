context("PlotCatalogToPdf.ID")

test_that("PlotCatalogToPdf.ID function is working properly", {
  catalog <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                       ref.genome = "GRCh37",
                       region = "genome", catalog.type = "counts")
  colnames(catalog) <- paste0("Biliary-AdenoCA", 1 : 35)
  out <- PlotCatalogToPdf(catalog, filename = "PlotCatID.test.pdf")
  expect_equal(out, TRUE)
  unlink("PlotCatID.test.pdf")
})

context("PlotCatalogToPdf.SNS96")

test_that("PlotCatalogToPdf.SNS96 function is working properly", {
  catalog <- ReadCatalog("testdata/regress.cat.sns.96.csv",
                         ref.genome = "GRCh37",
                         region = "genome", type = "counts")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog, filename = "PlotCatSNS96.test.pdf")
  expect_equal(out, TRUE)
  unlink("PlotCatSNS96.test.pdf")
  graphics.off()
})

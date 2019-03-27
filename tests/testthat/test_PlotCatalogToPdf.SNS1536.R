context("PlotCatalogToPdf.SNS1536")

test_that("PlotCatalogToPdf.SNS1536 function is working properly", {
  catalog <- ReadCatalog("testdata/regress.cat.sns.1536.csv",
                            ref.genome = "GRCh37",
                            region = "genome", catalog.type = "counts")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog, "PlotCatSNS1536.test.pdf")
  expect_equal(out, TRUE)
  unlink("PlotCatSNS1536.test.pdf")
  graphics.off()
})

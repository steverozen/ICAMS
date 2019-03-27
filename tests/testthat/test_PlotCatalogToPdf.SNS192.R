context("PlotCatalogToPdf.SNS192")

test_that("PlotCatalogToPdf.SNS192 function is working properly", {
  catalog <- ReadCatalog("testdata/regress.cat.sns.192.csv", ref.genome = "GRCh37",
                         region = "genome", catalog.type = "counts")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog, filename = "PlotCatSNS192.test.pdf")
  expect_equal(out, TRUE)
  unlink("PlotCatSNS192.test.pdf")
  graphics.off()
})

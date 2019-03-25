context("PlotCatalogToPdf.SNSClassStrandBias")

test_that("PlotCatalogToPdf.SNSClassStrandBias function is working properly", {
  catalog <- ReadCatalog("testdata/regress.cat.sns.192.csv", ref.genome = "GRCh37",
                         region = "genome", type = "counts")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog, filename = "PlotSNSClassStrandBias.test.pdf",
                     strandbias = TRUE)
  expect_equal(out, TRUE)
  unlink("PlotSNSClassStrandBias.test.pdf")
  graphics.off()
})

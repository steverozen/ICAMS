context("PlotCatalogToPdf.DNS136")

test_that("PlotCatalogToPdf.DNS136 function is working properly", {
  catalog <- ReadCatalog("testdata/regress.cat.dns.136.csv",
                           ref.genome = "GRCh37",
                           region = "genome", type = "counts")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog, filename = "PlotCatDNS136.test.pdf")
  expect_equal(out, TRUE)
  unlink("PlotCatDNS136.test.pdf")
  graphics.off()
})

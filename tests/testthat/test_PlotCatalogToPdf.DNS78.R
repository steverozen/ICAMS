context("PlotCatalogToPdf.DNS78")

test_that("PlotCatalogToPdf.DNS78 function is working properly", {
  catalog <- ReadCatalog("testdata/regress.cat.dns.78.csv",
                          ref.genome = "GRCh37",
                          region = "genome", catalog.type = "counts")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  out <- PlotCatalogToPdf(catalog, "PlotCatDNS78.test.pdf")
  expect_equal(out, TRUE)
  unlink("PlotCatDNS78.test.pdf")
  graphics.off()
})

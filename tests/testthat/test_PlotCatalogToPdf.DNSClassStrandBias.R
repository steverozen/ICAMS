context("PlotCatalogToPdf.DNSClassStrandBias")

test_that("PlotCatalogToPdf.DNSClassStrandBias function is working properly", {
  catalog <- ReadCatalog("testdata/regress.cat.dns.144.csv", ref.genome = "GRCh37",
                         region = "genome", type = "counts")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  out <- PlotCatalogToPdf(catalog,
                          filename = "PlotDNSClassStrandBias.test.pdf",
                          strandbias = TRUE)
  expect_equal(out, TRUE)
  unlink("PlotDNSClassStrandBias.test.pdf")
  graphics.off()
})

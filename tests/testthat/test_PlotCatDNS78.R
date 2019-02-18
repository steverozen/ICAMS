context("PlotCatDNS78")

test_that("PlotCatDNS78 function is working properly", {
  par(mar = c(6, 4, 6, 2))
  catalog <- ReadCatDNS78("testdata/regress.cat.dns.78.csv")
  out <- PlotCatDNS78(catalog[, 1, drop = FALSE], id = "test",
                      type = "density", abundance = abundance.2bp.GRCh37)
  expect_equal(out, TRUE)
})

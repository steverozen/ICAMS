context("PlotCatQUAD136")

test_that("PlotCatQUAD136 function is working properly", {
  par(oma = c(2, 2, 2, 0))
  catalog <- ReadCatDNS136("testdata/regress.cat.dns.136.csv")
  out <- PlotCatDNS136(catalog[, 1, drop = FALSE], id = "test",
                      type = "density", abundance = abundance.4bp.GRCh37)
  expect_equal(out, TRUE)
})

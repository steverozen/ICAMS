context("PlotDNSClassStrandBias")

test_that("PlotDNSClassStrandBias function is working properly", {
  par(mar = c(5, 8, 5, 1))
  catalog <- ReadCatDNS144("testdata/regress.cat.dns.144.csv")
  out <- PlotDNSClassStrandBias(catalog[, 1, drop = FALSE], "test")
  expect_equal(out, TRUE)
})

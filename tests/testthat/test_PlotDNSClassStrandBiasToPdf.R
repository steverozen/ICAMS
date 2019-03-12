context("PlotDNSClassStrandBiasToPdf")

test_that("PlotDNSClassStrandBiasToPdf function is working properly", {
  catalog <- ReadCatDNS144("testdata/regress.cat.dns.144.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat2 <- catalog[, 2, drop = FALSE]
  cat3 <- catalog[, 3, drop = FALSE]
  cat4 <- catalog[, 4, drop = FALSE]
  cat <- cbind(cat1, cat1, cat1, cat2, cat2, cat2, cat3, cat3, cat3,
               cat4, cat4, cat4)
  type <- rep(c("counts", "signature"), each = 3)
  out <- PlotDNSClassStrandBiasToPdf(cat, "PlotDNSClassStrandBias.test.pdf",
                                     type = rep(type, 2))
  expect_equal(out, TRUE)
  unlink("PlotDNSClassStrandBias.test.pdf")
  graphics.off()
})

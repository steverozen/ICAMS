context("PlotCatDNS78ToPdf")

test_that("PlotCatDNS78ToPdf function is working properly", {
  catalog <- ReadCatDNS78("testdata/regress.cat.dns.78.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat2 <- catalog[, 2, drop = FALSE]
  cat3 <- catalog[, 3, drop = FALSE]
  cat4 <- catalog[, 4, drop = FALSE]
  cat <- cbind(cat1, cat1, cat1, cat2, cat2, cat2, cat3, cat3, cat3,
               cat4, cat4, cat4)
  type <- c("counts", "signature", "density")
  out <-
    PlotCatDNS78ToPdf(cat, "PlotCatDNS78.test.pdf",
                      type = rep(type, 4),
                      abundance = abundance.2bp.genome.GRCh37)
  expect_equal(out, TRUE)
  unlink("PlotCatDNS78.test.pdf")
})

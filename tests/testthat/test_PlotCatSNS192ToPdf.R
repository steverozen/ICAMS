context("PlotCatSNS192ToPdf")

test_that("PlotCatSNS192ToPdf function is working properly", {
  catalog <- ReadCatSNS192("testdata/regress.cat.sns.192.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat2 <- catalog[, 2, drop = FALSE]
  cat3 <- catalog[, 3, drop = FALSE]
  cat4 <- catalog[, 4, drop = FALSE]
  cat <- cbind(cat1, cat1, cat2, cat2, cat3, cat3, cat4, cat4)
  type <- c("counts", "signature")
  out <- PlotCatSNS192ToPdf(cat, "PlotCatSNS192.test.pdf", type = rep(type, 4))
  expect_equal(out, TRUE)
})

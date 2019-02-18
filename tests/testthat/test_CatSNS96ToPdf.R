context("CatSNS96ToPdf")

test_that("CatSNS96ToPdf function is working properly", {
  catalog <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat2 <- catalog[, 2, drop = FALSE]
  cat3 <- catalog[, 3, drop = FALSE]
  cat4 <- catalog[, 4, drop = FALSE]
  cat <- cbind(cat1, cat1, cat1, cat2, cat2, cat2, cat3, cat3, cat3,
               cat4, cat4, cat4)
  type <- c("counts", "signature", "density")
  upper <- c(TRUE, rep(FALSE, 7), TRUE, rep(FALSE, 3))
  xlabels <- c(rep(FALSE, 7), TRUE, rep(FALSE, 3), TRUE)
  out <- CatSNS96ToPdf(cat, "PlotCatSNS96.test.pdf",
                    type = rep(type, 4), grid = FALSE, upper = upper,
                    xlabels = xlabels, abundance = .abundance.3bp.GRCh37)
  expect_equal(out, TRUE)
})

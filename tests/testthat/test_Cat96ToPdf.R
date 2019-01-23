context("Cat96ToPdf")

test_that("Cat96ToPdfNew function is working properly", {
  catalog <- ReadCat96("testdata/regress.cat.96.csv")
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
  out <- Cat96ToPdf(cat, "PlotCat96.test.pdf",
                    type = rep(type, 4), grid = FALSE, upper = upper,
                    xlabels = xlabels, abundance = .abundance.3bp)
  expect_equal(out, TRUE)
})

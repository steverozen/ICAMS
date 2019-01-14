context("Cat96ToPdfNew")

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
  out <- Cat96ToPdfNew(cat, "PlotCat96New.test.pdf",
                    type = rep(type, 4), abundance = .abundance.3bp)
  expect_equal(out, TRUE)
})

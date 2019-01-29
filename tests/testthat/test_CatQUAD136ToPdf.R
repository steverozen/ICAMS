context("CatQUAD136ToPdf")

test_that("CatQUAD136ToPdf function is working properly", {
  catalog <- ReadCatQUAD136("testdata/regress.cat.quad.136.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  type <- rep(c("density", "counts"), 2)
  out <- CatQUAD136ToPdf(catalog, "PlotCatQUAD136.test.pdf", type = type,
                         abundance = .abundance.4bp)
  expect_equal(out, TRUE)
})

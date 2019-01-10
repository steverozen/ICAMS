context("Cat1536toPdf")

test_that("Cat1536toPdf function is working properly", {
  catalog <- ReadCat1536("testdata/regress.cat.1536.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  out <-
    Cat1536ToPdf(catalog, "PlotCat1536.test.pdf", abundance = .abundance.5bp)
  expect_equal(out, TRUE)
})

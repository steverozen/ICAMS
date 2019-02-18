context("CatSNS1536toPdf")

test_that("CatSNS1536toPdf function is working properly", {
  catalog <- ReadCatSNS1536("testdata/regress.cat.sns.1536.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  out <-
    CatSNS1536ToPdf(catalog, "PlotCatSNS1536.test.pdf", abundance = .abundance.5bp.GRCh37)
  expect_equal(out, TRUE)
})

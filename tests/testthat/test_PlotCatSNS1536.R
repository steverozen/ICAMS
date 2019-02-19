context("PlotCatSNS1536")

test_that("PlotCatSNS1536 function is working properly", {
  catalog <- ReadCatSNS1536("testdata/regress.cat.sns.1536.csv")
  out <-
    PlotCatSNS1536(catalog[, 1, drop = FALSE],
                   "test",
                   abundance = abundance.5bp.genome.GRCh37)
  expect_equal(out, TRUE)
})

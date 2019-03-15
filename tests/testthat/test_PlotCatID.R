context("PlotCatID")

test_that("PlotCatID function is working properly", {
  par(mar = c(7, 4, 7, 3))
  catalog <- ReadCatID("testdata/BTSG_WGS_PCAWG.indels.csv")
  out <- PlotCatID(catalog[, 1, drop = FALSE], type = "counts", id = "test")
  expect_equal(out, TRUE)

  signature.catalog <- apply(catalog, MARGIN = 2, function (x) x / sum(x))
  out1 <- PlotCatID(signature.catalog[, 1, drop = FALSE],
                    type = "signature", id = "test")
  expect_equal(out1, TRUE)
})

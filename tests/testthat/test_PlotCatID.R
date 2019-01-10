context("PlotCatID")

test_that("PlotCatID function is working properly", {
  par(mar = c(7, 4, 7, 3))
  catalog <- ReadCatID("testdata/BTSG_WGS_PCAWG.indels.csv")
  out <- PlotCatID(catalog[, 1, drop = FALSE], id = "test")
  expect_equal(out, TRUE)
})

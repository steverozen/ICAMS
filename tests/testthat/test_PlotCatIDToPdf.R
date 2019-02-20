context("PlotCatIDToPdf")

test_that("PlotCatIDToPdf function is working properly", {
  catalog <- ReadCatID("testdata/BTSG_WGS_PCAWG.indels.csv")
  colnames(catalog) <- paste0("Biliary-AdenoCA", 1 : 35)
  cat1 <- catalog[, 1, drop = FALSE]
  cat2 <- catalog[, 2, drop = FALSE]
  cat3 <- catalog[, 3, drop = FALSE]
  cat4 <- catalog[, 4, drop = FALSE]
  cat <- cbind(cat1, cat1, cat2, cat2, cat3, cat3, cat4, cat4)
  type <- c("counts", "signature")
  out <- PlotCatIDToPdf(cat, "PlotCatID.test.pdf", type = rep(type, 4))
  expect_equal(out, TRUE)
})

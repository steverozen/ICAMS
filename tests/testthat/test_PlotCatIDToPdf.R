context("PlotCatIDToPdf")

test_that("PlotCatIDToPdf function is working properly", {
  catalog <- ReadCatID("testdata/BTSG_WGS_PCAWG.indels.csv")
  colnames(catalog) <- paste0("Biliary-AdenoCA", 1 : 35)

  cat1 <- catalog[, 1, drop = FALSE]
  cat1.sig <- apply(cat1, MARGIN = 2, function(x) x / sum(x))

  cat2 <- catalog[, 2, drop = FALSE]
  cat2.sig <- apply(cat2, MARGIN = 2, function(x) x / sum(x))

  cat3 <- catalog[, 3, drop = FALSE]
  cat3.sig <- apply(cat3, MARGIN = 2, function(x) x / sum(x))

  cat4 <- catalog[, 4, drop = FALSE]
  cat4.sig <- apply(cat4, MARGIN = 2, function(x) x / sum(x))

  cat <- cbind(cat1, cat1.sig,
               cat2, cat2.sig,
               cat3, cat3.sig,
               cat4, cat4.sig)

  type <- c("counts", "signature")
  out <- PlotCatIDToPdf(cat, name = "PlotCatID.test.pdf", type = rep(type, 4))
  expect_equal(out, TRUE)
  unlink("PlotCatID.test.pdf")
})

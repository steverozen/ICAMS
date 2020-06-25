context("Exposure related functions")

test_that("Read and write exposure", {
  file <- "testdata/tiny.exposure.csv"
  expect_warning(x <- ReadExposure(file))
  tfile <- tempfile()
  WriteExposure(x, tfile)
  reread.x <- ReadExposure(tfile)
  expect_equal(x, reread.x)
  
  x <- ReadExposure(file, check.names = FALSE)
  WriteExposure(x, tfile)
  x2 <- ReadExposure(tfile, check.names = FALSE)
  expect_equal(x, x2)
  
  file2 <- "testdata/tiny.exposure.dup.csv"
  expect_error(x <- ReadExposure(file2, check = FALSE))
})

test_that("PlotExposureInternal function", {
  exposure <- ReadExposure("testdata/synthetic.exposure.csv")
  PlotExposureInternal(SortExposure(exposure[, 1:30]), cex.legend = 0.5)
  PlotExposureInternal(SortExposure(exposure[ , 1:30]), plot.proportion = TRUE,
                       cex.legend = 0.5)
})
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

test_that("PlotExposure function", {
  exposure <- ReadExposure("testdata/synthetic.exposure.csv")
  
  old.par <- par(mfcol = c(2, 1), mar = c(2, 3, 3.5, 2), oma = c(2, 1, 0, 1))
  on.exit(par(old.par))
  PlotExposure(exposure = SortExposure(exposure[, 1:43]),
               main = "test", cex.main = 0.8, cex.legend = 0.35)
  
  par(old.par)
  PlotExposure(exposure = SortExposure(exposure[, 1:30]),  # Test a trick edge case
               main = "test1", cex.legend = 0.5)
  
  PlotExposure(exposure = SortExposure(exposure[, 3:6]), 
               samples.per.line = 4, cex.legend = 0.5)
  
  par(mfcol = c(2, 1), mar = c(2, 3, 1.5, 2), oma = c(2, 1, 0, 1))
  PlotExposure(exposure = SortExposure(exposure[, 1:43 ]),
               plot.proportion = TRUE, cex.legend = 0.3)
  
  par(old.par)
  PlotExposure(exposure = SortExposure(exposure[, 3:6]),
               samples.per.line = 4, plot.proportion = TRUE, 
               col = c("red", "white"), cex.legend = 0.5)
})

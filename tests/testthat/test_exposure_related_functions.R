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
  out <- PlotExposureInternal(SortExposure(exposure[, 1:30]), cex.legend = 0.5)
  expect_equal(out$plot.success, TRUE)
  
  out1 <- PlotExposureInternal(SortExposure(exposure[ , 1:30]), plot.proportion = TRUE,
                       cex.legend = 0.5)
  expect_equal(out1$plot.success, TRUE)
  
  out2 <- PlotExposureInternal(SortExposure(exposure[ , 3:6]), 
                               cex.legend = 0.5)
  expect_equal(out2$plot.success, TRUE)
})

test_that("PlotExposure function", {
  exposure <- ReadExposure("testdata/synthetic.exposure.csv")
  
  old.par <- par(mfcol = c(2, 1), mar = c(2, 3, 3.5, 2), oma = c(2, 1, 0, 1))
  on.exit(par(old.par))
  out <- PlotExposure(exposure = SortExposure(exposure[, 1:43]),
               main = "test", cex.main = 0.8, cex.legend = 0.35)
  expect_equal(out$plot.success, TRUE)
  
  par(old.par)
  out1 <- PlotExposure(exposure = SortExposure(exposure[, 1:30]),  # Test a trick edge case
               main = "test1", cex.legend = 0.5)
  expect_equal(out1$plot.success, TRUE)
  
  out2 <- PlotExposure(exposure = SortExposure(exposure[, 3:6]), 
               samples.per.line = 4, cex.legend = 0.5)
  expect_equal(out2$plot.success, TRUE)
  
  par(mfcol = c(2, 1), mar = c(2, 3, 1.5, 2), oma = c(2, 1, 0, 1))
  out3 <- PlotExposure(exposure = SortExposure(exposure[, 1:43 ]),
               plot.proportion = TRUE, cex.legend = 0.3)
  expect_equal(out3$plot.success, TRUE)
  
  par(old.par)
  out4 <- PlotExposure(exposure = SortExposure(exposure[, 3:6]),
               samples.per.line = 4, plot.proportion = TRUE, 
               col = c("red", "blue"), cex.legend = 0.5)
  expect_equal(out4$plot.success, TRUE)
})

test_that("PlotExposureToPdf function", {
  exposure <- ReadExposure("testdata/synthetic.exposure.csv")
  
  out <- PlotExposureToPdf(exposure,
                    file = file.path(tempdir(), "test.pdf"))
  expect_equal(out$plot.success, TRUE)
  
  out1 <- PlotExposureToPdf(exposure = SortExposure(exposure),
                    file = file.path(tempdir(), "test1.pdf"))
  expect_equal(out1$plot.success, TRUE)
  
  out2 <- PlotExposureToPdf(exposure = SortExposure(exposure[, 1:30]),
                    file = file.path(tempdir(), "test2.pdf"))  # Test a trick edge case
  expect_equal(out2$plot.success, TRUE)
  
  out3 <- PlotExposureToPdf(exposure = SortExposure(exposure[, 3:6]), 
                    file = file.path(tempdir(), "test3.pdf"), 
                    samples.per.line = 4)
  expect_equal(out3$plot.success, TRUE)
  
  out4 <- PlotExposureToPdf(exposure = SortExposure(exposure),
                    file = file.path(tempdir(), "test4.pdf"), 
                    plot.proportion = TRUE)
  expect_equal(out4$plot.success, TRUE)
  
  out5 <- PlotExposureToPdf(exposure = SortExposure(exposure[, 3:6]),
                    file = file.path(tempdir(), "test5.pdf"), 
                    samples.per.line = 4, plot.proportion = TRUE, 
                    col = c("red", "blue"))
  expect_equal(out5$plot.success, TRUE)
  unlink(file.path(tempdir(), paste0("test", c("", 1:5), ".pdf")))
  graphics.off()
})

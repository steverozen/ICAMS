context("PlotCatalog for matrix like object")

test_that("PlotCatalog for matrix like object", {
  mat96 <- matrix(data = rnbinom(96, mu = 100, size = 10), nrow = 96, ncol = 1)
  out96 <- PlotCatalog(mat96)
  expect_equal(out96$plot.success, TRUE)
  
  mat192 <- matrix(data = rnbinom(192, mu = 100, size = 10), 
                   nrow = 192, ncol = 1)
  out192 <- PlotCatalog(mat192)
  expect_equal(out192$plot.success, TRUE)
  
  mat1536 <- matrix(data = rnbinom(1536, mu = 100, size = 10), 
                    nrow = 1536, ncol = 1)
  out1536 <- PlotCatalog(mat1536)
  expect_equal(out1536$plot.success, TRUE)
  
  mat78 <- matrix(data = rnbinom(78, mu = 100, size = 10), 
                  nrow = 78, ncol = 1)
  out78 <- PlotCatalog(mat78)
  expect_equal(out78$plot.success, TRUE)
  
  mat136 <- round(matrix(data = rnbinom(136, mu = 100, size = 10), 
                         nrow = 136, ncol = 1))
  out136 <- PlotCatalog(mat136)
  expect_equal(out136$plot.success, TRUE)
  
  mat144 <- round(matrix(data = rnbinom(144, mu = 100, size = 10), 
                         nrow = 144, ncol = 1))
  out144 <- PlotCatalog(mat144)
  expect_equal(out144$plot.success, TRUE)
  
  mat83 <- round(matrix(data = rnbinom(83, mu = 100, size = 10), 
                        nrow = 83, ncol = 1))
  out83 <- PlotCatalog(mat83)
  expect_equal(out83$plot.success, TRUE)
  
})
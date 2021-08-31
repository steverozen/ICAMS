test_that("GeneratePlotPFMmatrix", {

  load("testdata/test_GeneratePlotPFMmatrix.Rdata")

  p1 <- tempfile(fileext = "pdf")
  
  test.INS.T.1.4.retval <- 
    GeneratePlotPFMmatrix(sequences = INS.T.1.4.sequences,
                          indel.class = "INS:T:1:4",
                          plot.dir=p1,
                          plot.title = "test.INS.T.1.4")

  p2 <- tempfile(fileext = "pdf")
  
  test.INS.T.1.0.retval <- 
    GeneratePlotPFMmatrix(sequences = INS.T.1.0.sequences,
                          indel.class = "INS:T:1:0",
                          plot.dir=p2,
                          plot.title = "test.INS.T.1.0")


  expect_equal(test.INS.T.1.4.retval, INS.T.1.4.retval)
  expect_equal(test.INS.T.1.0.retval, INS.T.1.0.retval)

})

test_that("errors", {
  expect_error(
    test.INS.T.1.4.retval <- GeneratePlotPFMmatrix(sequences = c(INS.T.1.4.sequences,INS.T.1.0.sequences),indel.class = "INS:T:1:4"))
})

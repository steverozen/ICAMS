test_that("GeneratePlotPFMmatrix", {

  load("testdata/test_GeneratePlotPFMmatrix.Rdata")

  test.INS.T.1.4.retval <- GeneratePlotPFMmatrix(sequences = INS.T.1.4.sequences,indel.class = "INS:T:1:4",plot.dir="test.INS.T.1.4.pdf",plot.title = "test.INS.T.1.4")

  test.INS.T.1.0.retval <- GeneratePlotPFMmatrix(sequences = INS.T.1.0.sequences,indel.class = "INS:T:1:0",plot.dir="test.INS.T.1.0.pdf",plot.title = "test.INS.T.1.0")


  expect_equal(test.INS.T.1.4.retval, INS.T.1.4.retval)
  expect_equal(test.INS.T.1.0.retval, INS.T.1.0.retval)

})

test_that("errors", {
  expect_error(
    test.INS.T.1.4.retval <- GeneratePlotPFMmatrix(sequences = c(INS.T.1.4.sequences,INS.T.1.0.sequences),indel.class = "INS:T:1:4"))
})
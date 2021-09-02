test_that("GeneratePlotPFMmatrix", {
  # Please see "./data-raw/code/create_testdata_for_ID_extended_sequences.R"
  # for generating the test data
  
  load("testdata/test_GeneratePlotPFMmatrix.Rdata")
  
  files <- file.path(tempdir(), paste0("test", 1:4, ".pdf"))
  
  test.INS.T.1.0.retval <- 
    GeneratePlotPFMmatrix(sequences = INS.T.1.0.sequences,
                          indel.class = "INS:T:1:0",
                          plot.dir = files[1],
                          plot.title = "test.INS.T.1.0")
  
  test.INS.T.1.4.retval <- 
    GeneratePlotPFMmatrix(sequences = INS.T.1.4.sequences,
                          indel.class = "INS:T:1:4",
                          plot.dir = files[2],
                          plot.title = "test.INS.T.1.4")

  test.DEL.C.1.0.retval <- 
    GeneratePlotPFMmatrix(sequences = DEL.C.1.0.sequences,
                          indel.class = "DEL:C:1:0",
                          plot.dir = files[3],
                          plot.title = "test.DEL.C.1.0")

  
  test.DEL.C.1.4.retval <- 
    GeneratePlotPFMmatrix(sequences = DEL.C.1.4.sequences,
                          indel.class = "DEL:C:1:4",
                          plot.dir = files[4],
                          plot.title = "test.DEL.C.1.4")
  
  expect_equal(test.INS.T.1.0.retval, INS.T.1.0.retval)
  expect_equal(test.INS.T.1.4.retval, INS.T.1.4.retval)
  expect_equal(test.DEL.C.1.0.retval, DEL.C.1.0.retval)
  expect_equal(test.DEL.C.1.4.retval, DEL.C.1.4.retval)
  
  sapply(files, FUN = unlink, recursive = TRUE)
})

test_that("errors", {
  expect_error(
    test.INS.T.1.4.retval <- 
      GeneratePlotPFMmatrix(sequences = c(INS.T.1.4.sequences, INS.T.1.0.sequences),
                            indel.class = "INS:T:1:4"))
})

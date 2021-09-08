test_that("HaplotypePlot", {
  # Please see "./data-raw/code/create_testdata_for_ID_extended_sequences.R"
  # for generating the test data

  load("testdata/test_HaplotypePlot.Rdata")

  test.INS.T.1.0.retval <-
    HaplotypePlot(sequences = INS.T.1.0.sequences,
                  indel.class = "INS:T:1:0",
                  title = "De novo insertion of 1T")

  test.INS.T.1.4.retval <-
    HaplotypePlot(sequences = INS.T.1.4.sequences,
                  indel.class = "INS:T:1:4",
                  title = "Insertion of 1T to 4Ts")

  test.DEL.C.1.0.retval <-
    HaplotypePlot(sequences = DEL.C.1.0.sequences,
                  indel.class = "DEL:C:1:0",
                  title = "Deletion of 1C from 1C")


  test.DEL.C.1.4.retval <-
    HaplotypePlot(sequences = DEL.C.1.4.sequences,
                  indel.class = "DEL:C:1:4",
                  title = "Deletion of 1C from 5Cs")

  expect_equal(test.INS.T.1.0.retval, INS.T.1.0.retval)
  expect_equal(test.INS.T.1.4.retval, INS.T.1.4.retval)
  expect_equal(test.DEL.C.1.0.retval, DEL.C.1.0.retval)
  expect_equal(test.DEL.C.1.4.retval, DEL.C.1.4.retval)

  plot(test.INS.T.1.0.retval)
  plot(test.INS.T.1.4.retval)
  plot(test.DEL.C.1.0.retval)
  plot(test.DEL.C.1.4.retval)
})

test_that("errors", {
  expect_error(
    test.INS.T.1.4.retval <-
      HaplotypePlot(sequences = c(INS.T.1.4.sequences, INS.T.1.0.sequences),
                    indel.class = "INS:T:1:4")
    )
})

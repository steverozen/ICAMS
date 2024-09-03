test_that("HaplotypePlot", {
  # Please see "./data-raw/code/create_testdata_for_ID_extended_sequences.R"
  # for generating the test data

  load("testdata/test_HaplotypePlot.Rdata")

  test.INS.T.1.0.retval <-
    HaplotypePlot(sequences = INS.T.1.0.sequences,
                  indel.class = "INS:T:1:0",
                  flank.length = 10,
                  title = "De novo insertion of 1T")

  test.INS.T.1.4.retval <-
    HaplotypePlot(sequences = INS.T.1.4.sequences,
                  indel.class = "INS:T:1:4",
                  flank.length = 6,
                  title = "Insertion of 1T to 4Ts")

  test.DEL.C.1.0.retval <-
    HaplotypePlot(sequences = DEL.C.1.0.sequences,
                  indel.class = "DEL:C:1:0",
                  flank.length = 10,
                  title = "Deletion of 1C from 1C")


  test.DEL.C.1.4.retval <-
    HaplotypePlot(sequences = DEL.C.1.4.sequences,
                  indel.class = "DEL:C:1:4",
                  flank.length = 6,
                  title = "Deletion of 1C from 5Cs")
  if (FALSE) {
    # Steve is not sure if it makes sense to compare
    # objects of class gg, ggplot. Possilby structure
    # of the objects change when ggplot2 changes.
    expect_equal(test.INS.T.1.0.retval, INS.T.1.0.retval)
    expect_equal(test.INS.T.1.4.retval, INS.T.1.4.retval)
    expect_equal(test.DEL.C.1.0.retval, DEL.C.1.0.retval)
    expect_equal(test.DEL.C.1.4.retval, DEL.C.1.4.retval)
  }
  
  expect_no_error(plot(test.INS.T.1.0.retval))
  expect_no_error(plot(test.INS.T.1.4.retval))
  expect_no_error(plot(test.DEL.C.1.0.retval))
  expect_no_error(plot(test.DEL.C.1.4.retval))
})

test_that("errors", {
  expect_error(
    test.INS.T.1.4.retval <-
      HaplotypePlot(sequences = c(INS.T.1.4.sequences, INS.T.1.0.sequences),
                    indel.class = "INS:T:1:4")
    )
})

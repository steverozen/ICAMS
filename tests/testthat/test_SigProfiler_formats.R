context("Reading and writing SigProfiler formats")

test_that("Read SigProfiler", {
  expect_warning(
  sp.cat <-
    ReadCatalog(
      "testdata/sigProfiler_ID_signatures.csv",
      strict = FALSE)
  )
  
})
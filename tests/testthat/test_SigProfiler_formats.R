context("Reading and writing SigProfiler formats")

test_that("Read SigProfiler", {
  expect_warning(
  sp.cat <-
    ReadCatalog(
      system.file(
        "tests/testthat/testdata/sigProfiler_ID_signatures.csv",
        package = "ICAMS",
        mustWork = TRUE),
      strict = FALSE)
  )
  
})
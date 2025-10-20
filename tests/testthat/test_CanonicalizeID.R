context("CanonicalizeID")
library(testthat)


test_that("CanonicalizeID", {
  source("c:/Users/steve/Documents/GitHub/ICAMS/inst/new_id_fns.R")
  dd = readr::read_csv(
    "c:/Users/steve/Documents/GitHub/ICAMS/tests/testthat/testdata/input_for_canonicalizeid.csv",
    col_types = "cccd"
  )

  xx = ICAMS:::CanonicalizeID(dd$context, dd$ref, dd$alt, dd$pos)
  expect_snapshot_value(xx, style = "deparse")
})

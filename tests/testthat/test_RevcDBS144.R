context("RevcDBS144")

test_that("RevcDBS144 function is working properly", {
  expect_equal(RevcDBS144(c("ACGT", "CTGA")), c("GTAC", "AGTC"))
})

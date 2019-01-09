context("RevcDNS144")

test_that("RevcDNS144 function is working properly", {
  expect_equal(RevcDNS144(c("ACGT", "CTGA")), c("GTAC", "AGTC"))
})

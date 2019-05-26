context("RevcSNS96")

test_that("RevcSNS96 function is working properly", {
  expect_equal(RevcSNS96(c("ATAG", "TGGA")), c("TATC", "CCAT"))
})

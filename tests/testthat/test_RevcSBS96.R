context("RevcSBS96")

test_that("RevcSBS96 function is working properly", {
  expect_equal(RevcSBS96(c("ATAG", "TGGA")), c("TATC", "CCAT"))
})

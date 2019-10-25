context("PyrPenta")

test_that("PyrPenta function", {
  to.test<- c("ACAAGA", "ACTAGA", "TTGCCA", "TTCGGA")
  expect_equal(PyrPenta(to.test), c("CTTGTT", "ACTAGA", "GGCAAT", "TTCGGA"))
})

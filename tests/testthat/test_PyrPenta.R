context("PyrPenta")

test_that("PyrPenta function is working properly", {
  to.test<- c("ACAAGA", "ACTAGA", "TTGCCA", "TTCGGA")
  expect_equal(PyrPenta(to.test), c("CTTGTT", "ACTAGA", "GGCAAT", "TTCGGA"))
})

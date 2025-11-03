test_that("justify_and_categorize_1_indel", {
  # the COSMIC-83 string should be "INS:C:1:0" (deletion of 1 C in a "repeat of 0 Cs")
  print(getwd())
  source("../../inst/xcategorize_1_justified_indel.R")
  x1 = ICAMS:::justify_and_categorize_1_indel(
    "AATCCC",
    "T",
    "TG",
    3
  )

  expect_snapshot(x1)
})

test_that("justify_and_categorize_1_indel", {
  # the COSMIC-83 string should be "INS:C:1:0" (deletion of 1 C in a "repeat of 0 Cs")

  x1 = ICAMS:::justify_and_categorize_1_indel(
    "AATCCC",
    "T",
    "TG",
    3
  )

  expect_snapshot(x1)
})

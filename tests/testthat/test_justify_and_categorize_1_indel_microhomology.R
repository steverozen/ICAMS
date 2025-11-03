test_that("justify_and_categorize_1_indel", {
  # the COSMIC-83 string should be "INS:C:1:0" (deletion of 1 C in a "repeat of 0 Cs")

  #  AAA[GC]GAA or
  #  AAAG[CG]AA
  #  AAA    GAA
  # insertion of CG
  # will justify to
  # AAA[GC]GAA
  # but orig repeat count will be 0

  source("../../inst/xcategorize_1_justified_indel.R")

  x1 = ICAMS:::justify_and_categorize_1_indel(
    "AAAGAACCC",
    "A",
    "AGC",
    3
  )

  expect_snapshot(x1)

  x2 = ICAMS:::justify_and_categorize_1_indel(
    "AAAGAACCC",
    "G",
    "GCG",
    4
  )

  expect_snapshot(x2)

  x3 = ICAMS:::justify_and_categorize_1_indel(
    "AAAGCGAACCC:",
    "GCG",
    "G",
    4
  )

  expect_snapshot(x3)

  x4 = ICAMS:::justify_and_categorize_1_indel(
    "AAAGCGAACCC:",
    "AGC",
    "A",
    3
  )

  expect_snapshot(x4)
})

test_that("justify_and_categorize_1_indel_insertion_R0", {
  # the COSMIC-83 string should be "INS:C:1:0" (deletion of 1 C in a "repeat of 0 Cs")

  #  AAA[GC]GAA or
  #  AAAG[CG]AA
  #  AAA    GAA
  # insertion of CG
  # will justify to
  # AAA[GC]GAA
  # but orig repeat count will be 0

  x1 = ICAMS:::justify_and_categorize_1_indel("CGCG", "G", "GAA", 3)
  expect_snapshot(x1)

  x2 = ICAMS:::justify_and_categorize_1_indel("CGCG", "G", "GAAA", 3)
  expect_snapshot(x2)

  x3 = ICAMS:::justify_and_categorize_1_indel("CGCG", "G", "GAAAA", 3)
  expect_snapshot(x3)

  x4 = ICAMS:::justify_and_categorize_1_indel("CGCG", "G", "GAT", 3)
  expect_snapshot(x4)

  x5 = ICAMS:::justify_and_categorize_1_indel("CGGG", "G", "GATC", 3) # Ins3:U3:R0
  expect_snapshot(x5)

  x6 = ICAMS:::justify_and_categorize_1_indel("CGGG", "G", "GATAT", 3)
  expect_snapshot(x6)

  x7 = ICAMS:::justify_and_categorize_1_indel("CGGG", "G", "GATTA", 3) # Ins4:U4:R0
  expect_snapshot(x7)

  x8 = ICAMS:::justify_and_categorize_1_indel("CGGG", "G", "GAAAAA", 3)
  expect_snapshot(x8)

  x9 = ICAMS:::justify_and_categorize_1_indel("CGAGG", "G", "GAA", 3)
  expect_snapshot(x9)

  x10 = ICAMS:::justify_and_categorize_1_indel("CGGG", "G", "GATATAT", 3)
  expect_snapshot(x10)

  x11 = ICAMS:::justify_and_categorize_1_indel("CGGG", "G", "GATTATT", 3)
  expect_snapshot(x11)
})

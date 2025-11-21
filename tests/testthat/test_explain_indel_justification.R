test_that("explain_indel_justification works with example cases", {
  # Test insertion case
  expect_message(
    explain_indel_justification("CAAAG", "CAAG", 3, 2, TRUE),
    "insertion"
  )

  # Test deletion cases
  expect_message(
    explain_indel_justification("ACTCTG", "ACTG", 3, 2, FALSE),
    "deletion"
  )

  expect_message(
    explain_indel_justification("ACTCTG", "ACTG", 4, 2, FALSE),
    "deletion"
  )

  expect_message(
    explain_indel_justification("CTCTG", "CTG", 3, 1, FALSE),
    "deletion"
  )
})

test_that("justify_indel examples work correctly", {
  # Example 1: deletion of "A" at position 2
  expect_snapshot(
    justify_indel("CAAAG", "CAAG", pos = 2, expected_delta = "A")
  )

  # Example 2: deletion of "A" at position 3
  expect_snapshot(
    justify_indel("CAAAG", "CAAG", pos = 3, expected_delta = "A")
  )

  # Example 3: deletion of "CA" at position 3, edge warning
  expect_snapshot(
    justify_indel("CACAG", "CAG", pos = 3, expected_delta = "CA")
  )

  # Deletion of "CA" at position 4, should move to the left
  expect_snapshot(
    justify_indel("TCACAG", "TCAG", pos = 4, expected_delta = "CA")
  )

  expect_snapshot(
    justify_indel("TCACAG", "TCAG", pos = 3, expected_delta = "AC")
  )

  expect_snapshot_value(
    justify_indel("TCACAG", "TCAG", pos = 3, expected_delta = "A")
  )
})

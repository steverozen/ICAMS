test_that("justify_many_indels works correctly", {
  # Read test data
  dd <- read.csv(
    "testdata/categorize_many_indels_test_input.csv"
  )

  # Store original POS values
  original_pos <- dd$POS

  # Run justify_many_indels
  result <- ICAMS:::justify_many_indels(dd, explain_indels = FALSE)

  # Check that pos_shift column was added
  expect_true("pos_shift" %in% colnames(result))

  # Check that pos_shift is always >= 0 (we only move left or stay)
  expect_true(all(result$pos_shift >= 0))

  # Check that POS was updated correctly (decremented by pos_shift)
  expect_equal(result$POS, original_pos - result$pos_shift)

  # Check that seq.context.width was updated correctly
  expect_equal(result$seq.context.width, dd$seq.context.width - result$pos_shift)

  # Check that categorization columns are present
  expect_true("COSMIC_83" %in% colnames(result))
  expect_true("Koh_89" %in% colnames(result))
  expect_true("Koh_476" %in% colnames(result))

  # Check that no rows were lost
  expect_equal(nrow(result), nrow(dd))
})

test_that("justify_many_indels handles empty input", {
  # Create empty dataframe with required columns
  empty_vcf <- data.frame(
    CHROM = character(0),
    POS = integer(0),
    REF = character(0),
    ALT = character(0),
    seq.context = character(0),
    seq.context.width = integer(0)
  )

  result <- ICAMS:::justify_many_indels(empty_vcf)

  expect_equal(nrow(result), 0)
  expect_true(is.data.frame(result))
})

test_that("justify_many_indels validates required columns", {
  # Create dataframe missing required column
  bad_vcf <- data.frame(
    CHROM = "1",
    POS = 100,
    REF = "AG",
    ALT = "A"
  )

  expect_error(
    ICAMS:::justify_many_indels(bad_vcf),
    "VCF missing required columns"
  )
})

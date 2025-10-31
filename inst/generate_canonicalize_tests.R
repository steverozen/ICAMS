# Temporary script to test Canonicalize1ID and capture results
# Run this in R console after loading the package with devtools::load_all()

library(testthat)
devtools::load_all()

# Test cases for Canonicalize1ID

cat("\n=== 1bp Deletions ===\n")

# 1bp deletion of C in a repeat (existing test)
result1 <- ICAMS:::Canonicalize1ID(
  "CTGTCCCCTTTCTTGAAGGTTTGGCTTCATATGCAATTCTAAGT",
  "G", "", 22
)
cat("1bp del C (repeat): ", result1, "\n")

# 1bp deletion of T with no repeat
result2 <- ICAMS:::Canonicalize1ID(
  "AAAAAATCCCCCC",
  "A", "", 7
)
cat("1bp del T (no repeat): ", result2, "\n")

# 1bp deletion of T in a repeat (TTT)
result3 <- ICAMS:::Canonicalize1ID(
  "AAAAAATTTCCCCCC",
  "A", "", 8
)
cat("1bp del T (in repeat): ", result3, "\n")

# 1bp deletion of A (should convert to T)
result4 <- ICAMS:::Canonicalize1ID(
  "CCCCAAAGGGGG",
  "T", "", 7
)
cat("1bp del A->T: ", result4, "\n")

# 1bp deletion of G (should convert to C)
result5 <- ICAMS:::Canonicalize1ID(
  "CCCCGGGAAAAA",
  "C", "", 7
)
cat("1bp del G->C: ", result5, "\n")

cat("\n=== Multi-bp Deletions in Repeats ===\n")

# 2bp deletion in repeat
result6 <- ICAMS:::Canonicalize1ID(
  "AAACACACACACGGG",
  "AC", "", 8
)
cat("2bp del in repeat: ", result6, "\n")

# 3bp deletion in repeat
result7 <- ICAMS:::Canonicalize1ID(
  "AAACAGCAGCAGCAGGGG",
  "CAG", "", 9
)
cat("3bp del in repeat: ", result7, "\n")

# 5bp deletion in repeat (should be 5+)
result8 <- ICAMS:::Canonicalize1ID(
  "AAAACTAGCACTAGCACTAGCACTAGCACTAGCAGGG",
  "CTAGC", "", 15
)
cat("5bp del in repeat: ", result8, "\n")

cat("\n=== Deletions with Microhomology ===\n")

# Using the example from FindDelMH documentation
# GAGAGG[CTAGAA]CTAGTT - has 4bp microhomology
result9 <- ICAMS:::Canonicalize1ID(
  "GGAGAGGCTAGAACTAGTTAAAAA",
  "CTAGAA", "", 7  # pos-1 for the deletion
)
cat("Deletion with MH (4bp): ", result9, "\n")

# 2bp deletion with no microhomology, no repeat
result10 <- ICAMS:::Canonicalize1ID(
  "AAAAAACGTTGGGGGG",
  "CG", "", 8
)
cat("2bp del no MH/repeat: ", result10, "\n")

cat("\n=== 1bp Insertions ===\n")

# 1bp insertion of T with no repeat
result11 <- ICAMS:::Canonicalize1ID(
  "AAAAAACCCCCC",
  "", "T", 6
)
cat("1bp ins T (no repeat): ", result11, "\n")

# 1bp insertion of T into existing T's
result12 <- ICAMS:::Canonicalize1ID(
  "AAAAAAATTTCCCCCC",
  "", "T", 7
)
cat("1bp ins T (in repeat): ", result12, "\n")

# 1bp insertion of A (should convert to T)
result13 <- ICAMS:::Canonicalize1ID(
  "CCCCCCGGGGGG",
  "", "A", 6
)
cat("1bp ins A->T: ", result13, "\n")

# 1bp insertion of G (should convert to C)
result14 <- ICAMS:::Canonicalize1ID(
  "TTTTTTAAAAAA",
  "", "G", 6
)
cat("1bp ins G->C: ", result14, "\n")

cat("\n=== Multi-bp Insertions ===\n")

# 2bp insertion in repeat context
result15 <- ICAMS:::Canonicalize1ID(
  "AAACACACGGG",
  "", "AC", 7
)
cat("2bp ins in repeat: ", result15, "\n")

# 3bp insertion in repeat
result16 <- ICAMS:::Canonicalize1ID(
  "AAACAGCAGCAGGGG",
  "", "CAG", 9
)
cat("3bp ins in repeat: ", result16, "\n")

# 2bp insertion no repeat
result17 <- ICAMS:::Canonicalize1ID(
  "AAAAAACCCCCC",
  "", "TG", 6
)
cat("2bp ins no repeat: ", result17, "\n")

# 5bp insertion (should be 5+)
result18 <- ICAMS:::Canonicalize1ID(
  "AAAAAACCCCCC",
  "", "TGCAT", 6
)
cat("5bp ins: ", result18, "\n")

cat("\n=== Edge Cases ===\n")

# Deletion with 5+ repeat count
result19 <- ICAMS:::Canonicalize1ID(
  "AAATTTTTTTTTCCCC",  # 9 T's total
  "A", "", 10
)
cat("1bp del with 5+ repeats: ", result19, "\n")

# Print all results in a format easy to copy to test
cat("\n=== Copy these to test file ===\n")
cat("test_that('Canonicalize1ID handles 1bp deletions', {\n")
cat("  expect_equal(ICAMS:::Canonicalize1ID('CTGTCCCCTTTCTTGAAGGTTTGGCTTCATATGCAATTCTAAGT', 'G', '', 22), '", result1, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAAAAATCCCCCC', 'A', '', 7), '", result2, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAAAAATTTCCCCCC', 'A', '', 8), '", result3, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('CCCCAAAGGGGG', 'T', '', 7), '", result4, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('CCCCGGGAAAAA', 'C', '', 7), '", result5, "')\n", sep="")
cat("})\n\n")

cat("test_that('Canonicalize1ID handles multi-bp deletions in repeats', {\n")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAACACACACACGGG', 'AC', '', 8), '", result6, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAACAGCAGCAGCAGGGG', 'CAG', '', 9), '", result7, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAAACTAGCACTAGCACTAGCACTAGCACTAGCAGGG', 'CTAGC', '', 15), '", result8, "')\n", sep="")
cat("})\n\n")

cat("test_that('Canonicalize1ID handles deletions with microhomology', {\n")
cat("  expect_equal(ICAMS:::Canonicalize1ID('GGAGAGGCTAGAACTAGTTAAAAA', 'CTAGAA', '', 7), '", result9, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAAAAACGTTGGGGGG', 'CG', '', 8), '", result10, "')\n", sep="")
cat("})\n\n")

cat("test_that('Canonicalize1ID handles 1bp insertions', {\n")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAAAAACCCCCC', '', 'T', 6), '", result11, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAAAAAATTTCCCCCC', '', 'T', 7), '", result12, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('CCCCCCGGGGGG', '', 'A', 6), '", result13, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('TTTTTTAAAAAA', '', 'G', 6), '", result14, "')\n", sep="")
cat("})\n\n")

cat("test_that('Canonicalize1ID handles multi-bp insertions', {\n")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAACACACGGG', '', 'AC', 7), '", result15, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAACAGCAGCAGGGG', '', 'CAG', 9), '", result16, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAAAAACCCCCC', '', 'TG', 6), '", result17, "')\n", sep="")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAAAAACCCCCC', '', 'TGCAT', 6), '", result18, "')\n", sep="")
cat("})\n\n")

cat("test_that('Canonicalize1ID handles edge cases', {\n")
cat("  expect_equal(ICAMS:::Canonicalize1ID('AAATTTTTTTTTCCCC', 'A', '', 10), '", result19, "')\n", sep="")
cat("})\n")

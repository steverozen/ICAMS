# Instructions for Creating Canonicalize1ID Tests

## Step 1: Run the test generation script

In R, from the ICAMS directory:

```r
source("generate_canonicalize_tests.R")
```

This will:
- Run multiple test cases for Canonicalize1ID
- Print the results with descriptions
- Generate formatted `expect_equal()` statements at the end

## Step 2: Copy the generated test code

The script will output formatted test code at the bottom. Copy these lines and paste them into the appropriate sections in `tests/testthat/test_Canonicalize1ID.R`.

## Step 3: Review and run the tests

```r
# Run just the Canonicalize1ID tests
testthat::test_file("tests/testthat/test_Canonicalize1ID.R")

# Or run all tests
devtools::test()
```

## What the tests cover

The generated tests include:

1. **1bp deletions**
   - Deletions with no repeats
   - Deletions in homopolymer repeats
   - Purine/pyrimidine conversion (A→T, G→C)

2. **Multi-bp deletions in repeats**
   - 2bp, 3bp deletions in tandem repeats
   - 5bp deletions (testing the "5+" classification)

3. **Deletions with microhomology**
   - Deletions with microhomology regions
   - Deletions without microhomology or repeats

4. **1bp insertions**
   - Insertions with no repeats
   - Insertions into existing repeats
   - Purine/pyrimidine conversion

5. **Multi-bp insertions**
   - Insertions in repeat contexts
   - Insertions with no repeat context
   - Large insertions (5+)

6. **Edge cases**
   - High repeat counts (5+)
   - Large indel sizes (5+)

## File locations

- **Test generation script**: `generate_canonicalize_tests.R`
- **Test file**: `tests/testthat/test_Canonicalize1ID.R`
- **Original draft tests**: `tests/testthat/xtest_Canonicalize1ID.R`

## Understanding the output format

The canonical representation format is:
- Deletions: `DEL:{base|repeats|MH}:{length}:{count}`
- Insertions: `INS:{base|repeats}:{length}:{count}`

Examples:
- `DEL:T:1:2` = 1bp deletion of T (converted from A) in 2 tandem repeats
- `DEL:repeats:3:1` = 3bp deletion in 1 repeat unit
- `DEL:MH:5:3` = 5bp deletion with 3bp microhomology
- `INS:C:1:0` = 1bp insertion of C with no repeats
- `INS:repeats:2:4` = 2bp insertion in 4 repeat units

## Notes

- The function automatically converts purines to pyrimidines (A→T, G→C)
- For deletions, `Canonicalize1ID` calls `Canonicalize1Del` with `pos + 1`
- For insertions, it calls `Canonicalize1INS` with the same `pos`
- Repeat counts of 5 or more are represented as "5+"
- Indel sizes of 5 or more are represented as "5+"

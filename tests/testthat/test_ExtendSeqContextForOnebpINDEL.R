
test_that("ExtendSeqContextForOnebpINDEL", {
  # Please see "./data-raw/code/create_testdata_for_ID_extended_sequences.R"
  # for generating the test data
  
  load("testdata/test_ExtendSeqContextForOnebpINDEL.Rdata")

  test.INS.T.1.0.sequences <- 
    SymmetricalContextsFor1BPIndel(test.indel.vcf, indel.class = "INS:T:1:0")
  test.INS.T.1.4.sequences <- 
    SymmetricalContextsFor1BPIndel(test.indel.vcf, indel.class = "INS:T:1:4")
  test.DEL.C.1.0.sequences <- 
    SymmetricalContextsFor1BPIndel(test.indel.vcf, indel.class = "DEL:C:1:0")
  test.DEL.C.1.4.sequences <- 
    SymmetricalContextsFor1BPIndel(test.indel.vcf, indel.class = "DEL:C:1:4")

  expect_equal(test.INS.T.1.0.sequences, INS.T.1.0.sequences)
  expect_equal(test.INS.T.1.4.sequences, INS.T.1.4.sequences)
  expect_equal(test.DEL.C.1.0.sequences, DEL.C.1.0.sequences)
  expect_equal(test.DEL.C.1.4.sequences, DEL.C.1.4.sequences)
})

test_that("errors", {
  expect_error(
    Get1BPIndelFlanks(sequence = "CTCTTTGTGAAGTAGCAACAGTTACATTTAAAATTTAAAACAC", 
                      ref = "T",  alt = "TT",  indel.class = "INS:T:1:2"))
# [1] "ACAGTTACATTT"
  expect_error(
    Get1BPIndelFlanks(sequence = "ATATAATTTGAATAAAAAGGTTTTGGAATACTGAAAAGTCTCC", 
                      ref = "T",  alt = "TT",  indel.class = "INS:T:1:4"))
# [1] "AGGTTTTGGAATAC"
  expect_error(
    Get1BPIndelFlanks(sequence = "CTGGCTAATTTTTTGTATTGTTTTTAGTAGAGACGGGGGTTTGC", 
                      ref = "TT",  alt = "T",  indel.class = "DEL:T:1:4"))
# [1] "TTGTTTTTAGTAGAG"
})

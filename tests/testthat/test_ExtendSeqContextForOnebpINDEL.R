
test_that("ExtendSeqContextForOnebpINDEL", {

  load("testdata/test_ExtendSeqContextForOnebpINDEL.Rdata")

  #This is to generate test.ins.vcf from a vcf file, but this does not need to run every time
  #test.ins.vcf <- ICAMS::ReadVCFs("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf")

  #test.ins.vcf <- ICAMS::VCFsToIDCatalogs(list(test.ins.vcf$Strelka.ID.GRCh37.s1),ref.genome="hg19",region = "genome",return.annotated.vcfs = T)

  #test.ins.vcf <- test.ins.vcf$annotated.vcfs[[1]]


  test.INS.T.1.0.sequences <- SymmetricalContextsFor1BPIndel(test.ins.vcf,indel.class = "INS:T:1:0")
  test.INS.T.1.4.sequences <- SymmetricalContextsFor1BPIndel(test.ins.vcf,indel.class = "INS:T:1:4")
  test.DEL.C.1.0.sequences <- SymmetricalContextsFor1BPIndel(test.ins.vcf,indel.class = "DEL:C:1:0")
  test.DEL.C.1.4.sequences <- SymmetricalContextsFor1BPIndel(test.ins.vcf,indel.class = "DEL:C:1:4")


  expect_equal(test.INS.T.1.0.sequences, INS.T.1.0.sequences)
  expect_equal(test.INS.T.1.4.sequences, INS.T.1.4.sequences)
  expect_equal(test.DEL.C.1.0.sequences, DEL.C.1.0.sequences)
  expect_equal(test.DEL.C.1.4.sequences, DEL.C.1.4.sequences)
})

test_that("errors", {
  expect_error(
    Get1BPIndelFlanks("CTCTTTGTGAAGTAGCAACAGTTACATTTAAAATTTAAAACAC", "T",  "TT",  "INS:T:1:2"))
# [1] "ACAGTTACATTT"
  expect_error(
    Get1BPIndelFlanks("ATATAATTTGAATAAAAAGGTTTTGGAATACTGAAAAGTCTCC", "T",  "TT",  "INS:T:1:4"))
# [1] "AGGTTTTGGAATAC"
  expect_error(
    Get1BPIndelFlanks("CTGGCTAATTTTTTGTATTGTTTTTAGTAGAGACGGGGGTTTGC", "TT",  "T",  "DEL:T:1:4"))
# [1] "TTGTTTTTAGTAGAG"
})

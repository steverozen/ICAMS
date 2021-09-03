# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

devtools::load_all()

test.indel.vcf <- 
  ReadVCFs("./tests/testthat/testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf")

test.indel.vcf <- 
  VCFsToIDCatalogs(test.indel.vcf,
                   ref.genome = "hg19",
                   region = "genome",
                   return.annotated.vcfs = TRUE)

test.indel.vcf <- test.indel.vcf$annotated.vcfs[[1]]

INS.T.1.0.sequences <- 
  SymmetricalContextsFor1BPIndel(test.indel.vcf, indel.class = "INS:T:1:0")
INS.T.1.4.sequences <- 
  SymmetricalContextsFor1BPIndel(test.indel.vcf, indel.class = "INS:T:1:4")
DEL.C.1.0.sequences <- 
  SymmetricalContextsFor1BPIndel(test.indel.vcf, indel.class = "DEL:C:1:0")
DEL.C.1.4.sequences <- 
  SymmetricalContextsFor1BPIndel(test.indel.vcf, indel.class = "DEL:C:1:4")

INS.T.1.0.retval <- 
  GeneratePlotPFMmatrix(sequences = INS.T.1.0.sequences,
                        indel.class = "INS:T:1:0")

INS.T.1.4.retval <- 
  GeneratePlotPFMmatrix(sequences = INS.T.1.4.sequences,
                        indel.class = "INS:T:1:4")

DEL.C.1.0.retval <-
  GeneratePlotPFMmatrix(sequences = DEL.C.1.0.sequences,
                        indel.class = "DEL:C:1:0")

DEL.C.1.4.retval <-
  GeneratePlotPFMmatrix(sequences = DEL.C.1.4.sequences,
                        indel.class = "DEL:C:1:4")

save(test.indel.vcf, INS.T.1.0.sequences, INS.T.1.4.sequences, 
     DEL.C.1.0.sequences, DEL.C.1.4.sequences, 
     file = "./tests/testthat/testdata/test_ExtendSeqContextForOnebpINDEL.Rdata")

save(INS.T.1.0.sequences, INS.T.1.4.sequences, DEL.C.1.0.sequences, 
     DEL.C.1.4.sequences, INS.T.1.0.retval, INS.T.1.4.retval, DEL.C.1.0.retval,
     DEL.C.1.4.retval,
     file = "./tests/testthat/testdata/test_GeneratePlotPFMmatrix.Rdata")
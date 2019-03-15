# Source this file from the top level directory.

cat(getwd(), "\n")

MakeTestDelVCF <- function() {
  return(
    data.frame(
      seq.context = c(
        "GAGGTATACATTGTGTTTACTTTTTCTATGTTTATGTACAATAGTAATATCTTTATAGTTATACTAACGTTATTAAAATAAGTAATTATATTAACTAAGTTTAGGACCAGTTTCTAGT",
        "GACCACTGAGAACCCAGGTTTTAGGCCCACCCCGGTACCAGGCCAGCCCCTGT",
        "AAGGTTTGGCTTCA",
        "ATTAAAATGGGGTT"),
      REF = c("ATAGTTATAC", "GCCCA", "TG", "AT"),
      ALT = c("A", "G", "T", "A"),
      seq.context.width = c(54, 24, 6, 6),
      stringsAsFactors = FALSE))
}

create.one.col.delete.test <-
  CreateOneColIDCatalog(MakeTestDelVCF(), NULL, trace = 2)
save(create.one.col.delete.test,
     file="tests/testthat/create_one_col_delete_test.Rdata")



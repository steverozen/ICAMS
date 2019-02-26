
context("Categorizing the mutation type of deletions")

test_that("FindMaxRepeatDel", {
  expect_equal(FindMaxRepeatDel("GATAGATAGATA", rep.unit.seq = "GATA", pos = 1), 2)
  expect_equal(FindMaxRepeatDel("GATAGATAGATA", rep.unit.seq = "GATA", pos = 5), 2)
  expect_equal(FindMaxRepeatDel("GATAGATAGATA", rep.unit.seq = "GATA", pos = 9), 2)
  expect_equal(FindMaxRepeatDel("CACAGATATATA", rep.unit.seq = "GATA", pos = 5), 0)
  expect_equal(FindMaxRepeatDel("CCCCCC", "C", pos = 3), 5)
  expect_error(FindMaxRepeatDel("CCCCCC", "T", pos = 3),
               "rep.unit.seq is not TRUE")
})

test_that("FindMaxRepeatIns", {
  tmp <- "xyaczt"
  expect_equal(FindMaxRepeatIns(tmp, "g", pos = 1), 0)
  expect_equal(FindMaxRepeatIns(tmp, "g", pos = 0), 0)
  expect_equal(FindMaxRepeatIns(tmp, "x", pos = 0), 1)
  expect_equal(FindMaxRepeatIns(tmp, "t", pos = 6), 1)
  expect_equal(FindMaxRepeatIns(tmp, "ac", pos = 2), 1)
  expect_equal(FindMaxRepeatIns(tmp, "ac", pos = 4), 1)
  for (i in c(1, 3, 5, 7, 9)) {
    expect_equal(FindMaxRepeatIns("gacacacacg", "ac", pos = i), 4)
  }
  for (i in c(0, 2, 4, 6, 8)) {
    expect_equal(FindMaxRepeatIns("acacacac", "ac", pos = i), 4)
  }
})

test_that("FindDelMH", {
  # GAGAGG[CTAGAA]CTAGTT
  #        ----   ----
  expect_equal(
    FindDelMH("GGAGAGGCTAGAACTAGTTAAAAA", "CTAGAA", 8, trace = 0),
    4)

  # GAGAGGC[TAGAAC]TAGTT
  #       * ---  * ---
  expect_equal(
    FindDelMH("GGAGAGGCTAGAACTAGTTAAAAA", "TAGAAC", 9, trace = 0),
    4)

  # AAGGCT[AGAACT]AGTTTT
  #     ** --  ** --
  expect_equal(
    FindDelMH("AAAGGCTAGAACTAGTTTTT", "AGAACT", 8, trace = 0),
    4)

  # GGCTA[GAACTA]GTT
  #   *** -  *** -
  expect_equal(
    FindDelMH("AAAGGCTAGAACTAGTTTTTT", "GAACTA", 9, trace = 0),
    4)

  # GGCTAG[AACTAG]TT
  #   ****   ****
  expect_equal(
    FindDelMH("AAAGGCTAGAACTAGTTTTTTT", "AACTAG", 10, trace = 0),
    4)

  # Cryptic repeat, return -1
  # TGACTA[GCTA]GTTAA
  #    *** -*** -
  expect_warning(
    FindDelMH("TGACTAGCTAGTTAA", "GCTA", 7, trace = 0),
    regexp = "unhandled cryptic repeat", fixed = TRUE)

  # Missed obvious repeat
  # AGATA[GATA]CCCCA
  #  **** ----
  expect_error(
    FindDelMH("AGATAGATACCCCA", "GATA", 6, trace = 0),
    "There is a repeated GATA to the left of the deleted GATA",
    fixed = TRUE)

  # Missed obvious repeat
  # ACCCCC[GATA]GATACCCCA
  #        **** ----
  expect_error(
    FindDelMH("ACCCCCGATAGATACCCCA", "GATA", 7, trace = 0),
    "There is a repeated GATA to the right of the deleted GATA",
    fixed = TRUE)

  # No microhomology at all
  # AAGATA[GATAG]CCCCAA
  #   **** ----
  expect_equal(
    FindDelMH("AAGATAGATAGCCCCAA", "GATAG", 7, trace = 0),
    0)

  # AAGATA[GGATA]CCCCAAA
  #   ****  ----
  expect_equal(
    FindDelMH("AAGATAGGATACCCCAAA", "GGATA", 7, trace = 0),
    4)
})

test_that("CreateOneColIDCatalog insertions", {
  MakeTestInsVCF <- function() {
    return(data.frame(
      seq.context=c("TTTTTTTTTTTTCGACCCCCCCCCCCC",
                    "TTTTTTTTTTTTGAACCCCCCCCCC",
                    "TTTTTTTTTTTTGCCCCCCCCCCCC",
                    "TTTTTTTTTTTTGCCCCCCCCCCCC"),
      REF=c("C", "G", "G", "G"),
      ALT=c("CGA", "GA", "GA", "GC"),
      seq.context.width=c(12, 12, 12, 12),
      stringsAsFactors = FALSE
    ))
  }
  load("create_one_col_insert_test.Rdata")
  expect_equal(
    ICAMS:::CreateOneColIDCatalog(MakeTestInsVCF(), NULL, trace = 0),
    create.one.col.insert.test)
})

test_that("CreateOneColIDCatalog deletions", {
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
  load("create_one_col_delete_test.Rdata")
  expect_equal(
    ICAMS:::CreateOneColIDCatalog(MakeTestDelVCF(), NULL, trace = 0),
    create.one.col.delete.test)
})


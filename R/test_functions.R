#' This function is to test whether the predefined functions
#' are working correctly to produce the desired SNS catalogs from Strelka VCF.
#'
#' @export
TestStrelkaSNSCatalog <- function() {
  expected.cat96 <-
    structure(
      c(25L, 77L, 1L, 38L, 59L, 28L, 8L, 71L, 46L, 191L,
        6L, 28L, 26L, 40L, 3L, 67L, 15L, 10L, 2L, 19L, 7L, 6L, 4L, 23L,
        7L, 82L, 5L, 25L, 12L, 11L, 2L, 71L, 26L, 15L, 15L, 14L, 103L,
        363L, 36L, 275L, 42L, 49L, 21L, 25L, 18L, 53L, 8L, 52L, 15L,
        13L, 8L, 23L, 100L, 82L, 138L, 109L, 7L, 10L, 13L, 10L, 12L,
        11L, 6L, 12L, 37L, 16L, 24L, 20L, 23L, 51L, 30L, 47L, 20L, 53L,
        13L, 9L, 18L, 14L, 9L, 14L, 10L, 8L, 27L, 8L, 8L, 11L, 18L, 5L,
        4L, 18L, 67L, 6L, 10L, 2L, 25L, 12L), .Dim = c(96L, 1L),
      .Dimnames =
        list(
          c("ACAA", "ACCA", "ACGA", "ACTA", "CCAA", "CCCA", "CCGA",
            "CCTA", "GCAA", "GCCA", "GCGA", "GCTA", "TCAA", "TCCA", "TCGA",
            "TCTA", "ACAG", "ACCG", "ACGG", "ACTG", "CCAG", "CCCG", "CCGG",
            "CCTG", "GCAG", "GCCG", "GCGG", "GCTG", "TCAG", "TCCG", "TCGG",
            "TCTG", "ACAT", "ACCT", "ACGT", "ACTT", "CCAT", "CCCT", "CCGT",
            "CCTT", "GCAT", "GCCT", "GCGT", "GCTT", "TCAT", "TCCT", "TCGT",
            "TCTT", "ATAA", "ATCA", "ATGA", "ATTA", "CTAA", "CTCA", "CTGA",
            "CTTA", "GTAA", "GTCA", "GTGA", "GTTA", "TTAA", "TTCA", "TTGA",
            "TTTA", "ATAC", "ATCC", "ATGC", "ATTC", "CTAC", "CTCC", "CTGC",
            "CTTC", "GTAC", "GTCC", "GTGC", "GTTC", "TTAC", "TTCC", "TTGC",
            "TTTC", "ATAG", "ATCG", "ATGG", "ATTG", "CTAG", "CTCG", "CTGG",
            "CTTG", "GTAG", "GTCG", "GTGG", "GTTG", "TTAG", "TTCG", "TTGG",
            "TTTG"), NULL))

  expected.cat192 <-
    structure(
      c(0L, 0L, 1L, 2L, 3L, 14L, 2L, 4L, 0L, 7L, 7L, 2L,
        0L, 1L, 2L, 3L, 5L, 2L, 18L, 3L, 1L, 2L, 10L, 8L, 2L, 16L, 19L,
        1L, 4L, 6L, 6L, 7L, 3L, 1L, 38L, 4L, 2L, 1L, 50L, 2L, 1L, 2L,
        19L, 2L, 3L, 1L, 26L, 6L, 3L, 15L, 0L, 6L, 10L, 4L, 2L, 14L,
        12L, 43L, 1L, 5L, 4L, 4L, 2L, 13L, 6L, 0L, 0L, 6L, 1L, 1L, 0L,
        2L, 0L, 15L, 2L, 7L, 2L, 1L, 1L, 11L, 5L, 1L, 6L, 1L, 17L, 57L,
        10L, 42L, 9L, 8L, 4L, 3L, 2L, 13L, 1L, 13L, 12L, 7L, 92L, 5L,
        0L, 2L, 10L, 3L, 17L, 11L, 111L, 3L, 6L, 11L, 30L, 8L, 14L, 3L,
        8L, 2L, 1L, 0L, 1L, 1L, 4L, 13L, 1L, 1L, 2L, 2L, 0L, 4L, 15L,
        3L, 17L, 11L, 0L, 1L, 3L, 0L, 11L, 53L, 9L, 22L, 8L, 9L, 15L,
        6L, 4L, 2L, 3L, 4L, 18L, 12L, 27L, 18L, 2L, 2L, 1L, 3L, 1L, 5L,
        2L, 4L, 11L, 5L, 3L, 6L, 6L, 10L, 7L, 7L, 5L, 6L, 3L, 2L, 2L,
        3L, 2L, 1L, 2L, 1L, 10L, 1L, 1L, 1L, 7L, 1L, 1L, 2L, 21L, 2L,
        2L, 1L, 5L, 3L),
      .Dim = c(192L, 1L),
      .Dimnames =
        list(
          c("AAAC",
            "AACC", "AAGC", "AATC", "CAAC", "CACC", "CAGC", "CATC", "GAAC",
            "GACC", "GAGC", "GATC", "TAAC", "TACC", "TAGC", "TATC", "AAAG",
            "AACG", "AAGG", "AATG", "CAAG", "CACG", "CAGG", "CATG", "GAAG",
            "GACG", "GAGG", "GATG", "TAAG", "TACG", "TAGG", "TATG", "AAAT",
            "AACT", "AAGT", "AATT", "CAAT", "CACT", "CAGT", "CATT", "GAAT",
            "GACT", "GAGT", "GATT", "TAAT", "TACT", "TAGT", "TATT", "ACAA",
            "ACCA", "ACGA", "ACTA", "CCAA", "CCCA", "CCGA", "CCTA", "GCAA",
            "GCCA", "GCGA", "GCTA", "TCAA", "TCCA", "TCGA", "TCTA", "ACAG",
            "ACCG", "ACGG", "ACTG", "CCAG", "CCCG", "CCGG", "CCTG", "GCAG",
            "GCCG", "GCGG", "GCTG", "TCAG", "TCCG", "TCGG", "TCTG", "ACAT",
            "ACCT", "ACGT", "ACTT", "CCAT", "CCCT", "CCGT", "CCTT", "GCAT",
            "GCCT", "GCGT", "GCTT", "TCAT", "TCCT", "TCGT", "TCTT", "AGAA",
            "AGCA", "AGGA", "AGTA", "CGAA", "CGCA", "CGGA", "CGTA", "GGAA",
            "GGCA", "GGGA", "GGTA", "TGAA", "TGCA", "TGGA", "TGTA", "AGAC",
            "AGCC", "AGGC", "AGTC", "CGAC", "CGCC", "CGGC", "CGTC", "GGAC",
            "GGCC", "GGGC", "GGTC", "TGAC", "TGCC", "TGGC", "TGTC", "AGAT",
            "AGCT", "AGGT", "AGTT", "CGAT", "CGCT", "CGGT", "CGTT", "GGAT",
            "GGCT", "GGGT", "GGTT", "TGAT", "TGCT", "TGGT", "TGTT", "ATAA",
            "ATCA", "ATGA", "ATTA", "CTAA", "CTCA", "CTGA", "CTTA", "GTAA",
            "GTCA", "GTGA", "GTTA", "TTAA", "TTCA", "TTGA", "TTTA", "ATAC",
            "ATCC", "ATGC", "ATTC", "CTAC", "CTCC", "CTGC", "CTTC", "GTAC",
            "GTCC", "GTGC", "GTTC", "TTAC", "TTCC", "TTGC", "TTTC", "ATAG",
            "ATCG", "ATGG", "ATTG", "CTAG", "CTCG", "CTGG", "CTTG", "GTAG",
            "GTCG", "GTGG", "GTTG", "TTAG", "TTCG", "TTGG", "TTTG"), NULL))

  expected.cat1536 <-
    structure(
      c(5L, 0L, 4L, 1L, 6L, 10L, 0L, 12L, 0L, 0L, 0L, 1L,
        5L, 0L, 5L, 1L, 2L, 7L, 2L, 5L, 2L, 2L, 3L, 4L, 2L, 1L, 0L, 0L,
        5L, 4L, 6L, 6L, 8L, 3L, 8L, 3L, 22L, 22L, 1L, 46L, 1L, 0L, 1L,
        0L, 1L, 1L, 2L, 5L, 0L, 0L, 3L, 0L, 3L, 0L, 1L, 4L, 1L, 0L, 0L,
        0L, 4L, 3L, 6L, 9L, 0L, 2L, 2L, 0L, 4L, 9L, 0L, 12L, 0L, 0L,
        0L, 0L, 4L, 2L, 5L, 3L, 7L, 4L, 5L, 5L, 2L, 0L, 1L, 4L, 0L, 0L,
        2L, 0L, 2L, 3L, 7L, 11L, 0L, 0L, 2L, 0L, 3L, 5L, 0L, 11L, 1L,
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 4L, 3L, 1L, 1L, 7L, 1L, 7L, 0L,
        0L, 1L, 0L, 7L, 2L, 4L, 6L, 4L, 0L, 2L, 3L, 1L, 1L, 1L, 0L, 0L,
        0L, 0L, 0L, 1L, 1L, 2L, 0L, 5L, 2L, 1L, 1L, 1L, 1L, 0L, 0L, 0L,
        0L, 0L, 1L, 1L, 1L, 2L, 2L, 3L, 3L, 5L, 2L, 2L, 5L, 2L, 10L,
        1L, 0L, 0L, 2L, 1L, 2L, 8L, 0L, 1L, 1L, 0L, 1L, 5L, 1L, 0L, 2L,
        0L, 1L, 0L, 0L, 2L, 4L, 2L, 2L, 0L, 0L, 2L, 0L, 5L, 7L, 0L, 9L,
        0L, 0L, 0L, 0L, 1L, 2L, 4L, 2L, 1L, 4L, 0L, 8L, 3L, 3L, 1L, 1L,
        1L, 1L, 0L, 0L, 5L, 1L, 5L, 10L, 4L, 0L, 2L, 3L, 10L, 23L, 0L,
        29L, 0L, 0L, 0L, 0L, 2L, 0L, 2L, 4L, 2L, 1L, 5L, 1L, 1L, 0L,
        1L, 6L, 0L, 0L, 0L, 0L, 3L, 3L, 5L, 5L, 1L, 0L, 1L, 0L, 0L, 0L,
        0L, 2L, 1L, 0L, 1L, 0L, 0L, 1L, 3L, 2L, 0L, 1L, 0L, 0L, 0L, 1L,
        0L, 0L, 0L, 1L, 0L, 0L, 1L, 2L, 1L, 1L, 0L, 1L, 2L, 2L, 19L,
        17L, 3L, 26L, 0L, 0L, 2L, 0L, 1L, 5L, 2L, 5L, 1L, 0L, 1L, 0L,
        1L, 1L, 0L, 2L, 0L, 0L, 0L, 0L, 3L, 5L, 14L, 5L, 1L, 0L, 2L,
        1L, 0L, 2L, 0L, 1L, 0L, 0L, 0L, 0L, 3L, 0L, 2L, 1L, 1L, 0L, 1L,
        1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 2L, 1L, 1L, 0L, 0L, 0L,
        0L, 1L, 0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 1L, 1L,
        1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 4L, 3L, 3L, 5L, 3L, 2L, 2L,
        1L, 2L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 1L,
        0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 4L, 0L, 0L, 0L, 0L,
        2L, 1L, 2L, 0L, 6L, 1L, 0L, 0L, 1L, 0L, 3L, 3L, 2L, 1L, 0L, 0L,
        0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 9L, 3L, 0L, 0L, 0L,
        1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 0L,
        0L, 2L, 1L, 0L, 1L, 0L, 0L, 2L, 0L, 4L, 0L, 3L, 1L, 0L, 0L, 0L,
        0L, 1L, 1L, 0L, 2L, 0L, 1L, 0L, 0L, 3L, 1L, 0L, 0L, 3L, 0L, 1L,
        0L, 0L, 2L, 0L, 3L, 0L, 0L, 0L, 0L, 4L, 2L, 6L, 4L, 3L, 0L, 2L,
        3L, 2L, 1L, 0L, 2L, 1L, 0L, 2L, 1L, 2L, 0L, 2L, 0L, 1L, 2L, 3L,
        1L, 4L, 0L, 2L, 7L, 1L, 0L, 1L, 4L, 4L, 9L, 4L, 9L, 5L, 3L, 2L,
        5L, 5L, 1L, 2L, 10L, 0L, 0L, 1L, 0L, 3L, 2L, 1L, 3L, 1L, 1L,
        1L, 1L, 0L, 4L, 0L, 5L, 1L, 0L, 0L, 0L, 2L, 2L, 3L, 2L, 2L, 0L,
        2L, 3L, 0L, 2L, 0L, 1L, 2L, 0L, 0L, 0L, 1L, 2L, 1L, 0L, 6L, 13L,
        12L, 12L, 55L, 0L, 16L, 35L, 0L, 2L, 3L, 6L, 6L, 16L, 11L, 22L,
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 2L, 1L, 0L, 0L, 1L, 0L,
        0L, 2L, 0L, 2L, 2L, 9L, 0L, 11L, 3L, 1L, 1L, 0L, 1L, 4L, 4L,
        7L, 1L, 1L, 0L, 2L, 1L, 2L, 0L, 1L, 0L, 2L, 2L, 0L, 0L, 1L, 0L,
        1L, 0L, 0L, 0L, 0L, 2L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 3L,
        0L, 1L, 1L, 3L, 2L, 3L, 1L, 0L, 5L, 0L, 1L, 5L, 2L, 1L, 2L, 1L,
        2L, 0L, 1L, 2L, 3L, 0L, 2L, 0L, 2L, 0L, 0L, 0L, 0L, 1L, 2L, 2L,
        5L, 1L, 2L, 2L, 2L, 1L, 1L, 0L, 1L, 1L, 1L, 0L, 3L, 0L, 1L, 1L,
        2L, 10L, 11L, 15L, 17L, 79L, 52L, 10L, 96L, 3L, 5L, 2L, 4L, 31L,
        56L, 63L, 38L, 6L, 5L, 5L, 4L, 1L, 8L, 0L, 11L, 2L, 2L, 2L, 3L,
        2L, 4L, 0L, 3L, 0L, 0L, 2L, 2L, 4L, 7L, 0L, 7L, 1L, 0L, 0L, 1L,
        3L, 1L, 1L, 12L, 0L, 1L, 2L, 1L, 0L, 0L, 2L, 8L, 3L, 1L, 0L,
        2L, 3L, 0L, 1L, 2L, 12L, 15L, 25L, 21L, 8L, 6L, 2L, 9L, 20L,
        16L, 26L, 19L, 19L, 8L, 22L, 13L, 2L, 0L, 1L, 0L, 1L, 2L, 0L,
        2L, 3L, 0L, 1L, 2L, 0L, 1L, 0L, 0L, 0L, 2L, 2L, 0L, 1L, 0L, 0L,
        0L, 0L, 1L, 0L, 2L, 1L, 0L, 1L, 1L, 1L, 0L, 1L, 2L, 0L, 1L, 0L,
        2L, 1L, 0L, 0L, 0L, 2L, 1L, 0L, 2L, 3L, 2L, 6L, 1L, 10L, 9L,
        2L, 11L, 5L, 3L, 13L, 6L, 8L, 6L, 3L, 5L, 1L, 0L, 0L, 0L, 0L,
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 2L, 0L,
        0L, 0L, 2L, 0L, 0L, 0L, 0L, 1L, 2L, 1L, 2L, 2L, 0L, 1L, 0L, 0L,
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 3L, 0L, 1L, 1L, 3L, 2L, 3L, 1L,
        1L, 0L, 3L, 1L, 1L, 3L, 3L, 1L, 2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L,
        0L, 0L, 0L, 0L, 2L, 2L, 1L, 1L, 1L, 0L, 2L, 0L, 2L, 0L, 1L, 1L,
        0L, 0L, 1L, 0L, 0L, 2L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 2L, 1L, 0L,
        0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 4L,
        3L, 2L, 11L, 6L, 5L, 8L, 3L, 5L, 4L, 6L, 5L, 1L, 0L, 2L, 0L,
        1L, 0L, 0L, 4L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 2L, 1L, 0L, 0L, 1L,
        3L, 2L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 6L, 2L, 2L, 3L,
        4L, 1L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 0L, 1L, 4L, 0L, 1L, 3L, 1L,
        3L, 4L, 0L, 1L, 3L, 0L, 0L, 0L, 5L, 1L, 2L, 2L, 0L, 1L, 1L, 2L,
        13L, 2L, 1L, 3L, 2L, 1L, 2L, 0L, 0L, 0L, 3L, 0L, 1L, 2L, 0L,
        0L, 3L, 0L, 0L, 3L, 0L, 1L, 0L, 0L, 0L, 0L, 2L, 0L, 2L, 2L, 4L,
        2L, 2L, 3L, 0L, 1L, 0L, 0L, 4L, 4L, 0L, 1L, 0L, 2L, 3L, 0L, 1L,
        1L, 7L, 19L, 1L, 3L, 3L, 6L, 3L, 3L, 1L, 1L, 1L, 5L, 0L, 0L,
        0L, 1L, 2L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 2L, 1L,
        1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 2L, 0L, 1L, 1L, 3L, 1L, 1L, 0L,
        4L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 4L, 1L, 0L, 0L, 0L, 0L, 1L, 1L,
        1L, 0L, 0L, 1L, 0L, 0L, 3L, 0L, 0L, 0L, 2L, 3L, 1L, 1L, 2L, 0L,
        0L, 1L, 12L, 0L, 3L, 4L, 1L, 3L, 1L, 1L, 0L, 0L, 0L, 0L, 0L,
        0L, 2L, 1L, 0L, 4L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 3L,
        1L, 2L, 2L, 1L, 0L, 0L, 2L, 3L, 2L, 1L, 0L, 2L, 2L, 1L, 5L, 1L,
        3L, 5L, 1L, 1L, 5L, 0L, 6L, 0L, 2L, 5L, 2L, 7L, 4L, 3L, 8L, 5L,
        0L, 4L, 3L, 6L, 1L, 1L, 4L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 2L, 3L,
        4L, 0L, 1L, 0L, 1L, 0L, 2L, 4L, 1L, 0L, 1L, 1L, 3L, 1L, 0L, 1L,
        1L, 1L, 1L, 0L, 0L, 0L, 2L, 3L, 2L, 3L, 2L, 2L, 1L, 0L, 0L, 1L,
        1L, 0L, 1L, 0L, 1L, 1L, 0L, 4L, 1L, 3L, 2L, 2L, 0L, 0L, 1L, 0L,
        0L, 0L, 0L, 1L, 0L, 0L, 5L, 2L, 1L, 3L, 0L, 0L, 1L, 3L, 0L, 3L,
        0L, 0L, 1L, 0L, 0L, 0L, 0L, 3L, 2L, 2L, 0L, 3L, 1L, 0L, 3L, 1L,
        0L, 0L, 1L, 0L, 1L, 0L, 1L, 3L, 2L, 4L, 1L, 0L, 0L, 0L, 2L, 0L,
        0L, 0L, 0L, 0L, 2L, 1L, 2L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 1L, 1L,
        0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L,
        0L, 2L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 3L, 1L, 1L, 0L, 0L, 0L, 0L,
        0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
        0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 0L, 0L, 0L, 0L, 0L, 0L,
        1L, 1L, 0L, 3L, 1L, 0L, 2L, 4L, 0L, 52L, 0L, 0L, 1L, 0L, 0L,
        0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 4L, 0L, 1L, 1L, 1L, 0L,
        2L, 0L, 1L, 0L, 2L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 1L,
        1L, 1L, 1L, 1L, 2L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L,
        0L, 0L, 1L, 0L, 2L, 1L, 0L, 1L, 0L, 0L, 3L, 1L, 0L, 0L, 0L, 1L,
        1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 2L, 1L, 2L, 1L, 0L, 0L, 0L
      ),
      .Dim = c(1536L, 1L),
      .Dimnames =
        list(c("AACAAA", "AACACA",
               "AACAGA", "AACATA", "AACCAA", "AACCCA", "AACCGA", "AACCTA", "AACGAA",
               "AACGCA", "AACGGA", "AACGTA", "AACTAA", "AACTCA", "AACTGA", "AACTTA",
               "ACCAAA", "ACCACA", "ACCAGA", "ACCATA", "ACCCAA", "ACCCCA", "ACCCGA",
               "ACCCTA", "ACCGAA", "ACCGCA", "ACCGGA", "ACCGTA", "ACCTAA", "ACCTCA",
               "ACCTGA", "ACCTTA", "AGCAAA", "AGCACA", "AGCAGA", "AGCATA", "AGCCAA",
               "AGCCCA", "AGCCGA", "AGCCTA", "AGCGAA", "AGCGCA", "AGCGGA", "AGCGTA",
               "AGCTAA", "AGCTCA", "AGCTGA", "AGCTTA", "ATCAAA", "ATCACA", "ATCAGA",
               "ATCATA", "ATCCAA", "ATCCCA", "ATCCGA", "ATCCTA", "ATCGAA", "ATCGCA",
               "ATCGGA", "ATCGTA", "ATCTAA", "ATCTCA", "ATCTGA", "ATCTTA", "CACAAA",
               "CACACA", "CACAGA", "CACATA", "CACCAA", "CACCCA", "CACCGA", "CACCTA",
               "CACGAA", "CACGCA", "CACGGA", "CACGTA", "CACTAA", "CACTCA", "CACTGA",
               "CACTTA", "CCCAAA", "CCCACA", "CCCAGA", "CCCATA", "CCCCAA", "CCCCCA",
               "CCCCGA", "CCCCTA", "CCCGAA", "CCCGCA", "CCCGGA", "CCCGTA", "CCCTAA",
               "CCCTCA", "CCCTGA", "CCCTTA", "CGCAAA", "CGCACA", "CGCAGA", "CGCATA",
               "CGCCAA", "CGCCCA", "CGCCGA", "CGCCTA", "CGCGAA", "CGCGCA", "CGCGGA",
               "CGCGTA", "CGCTAA", "CGCTCA", "CGCTGA", "CGCTTA", "CTCAAA", "CTCACA",
               "CTCAGA", "CTCATA", "CTCCAA", "CTCCCA", "CTCCGA", "CTCCTA", "CTCGAA",
               "CTCGCA", "CTCGGA", "CTCGTA", "CTCTAA", "CTCTCA", "CTCTGA", "CTCTTA",
               "GACAAA", "GACACA", "GACAGA", "GACATA", "GACCAA", "GACCCA", "GACCGA",
               "GACCTA", "GACGAA", "GACGCA", "GACGGA", "GACGTA", "GACTAA", "GACTCA",
               "GACTGA", "GACTTA", "GCCAAA", "GCCACA", "GCCAGA", "GCCATA", "GCCCAA",
               "GCCCCA", "GCCCGA", "GCCCTA", "GCCGAA", "GCCGCA", "GCCGGA", "GCCGTA",
               "GCCTAA", "GCCTCA", "GCCTGA", "GCCTTA", "GGCAAA", "GGCACA", "GGCAGA",
               "GGCATA", "GGCCAA", "GGCCCA", "GGCCGA", "GGCCTA", "GGCGAA", "GGCGCA",
               "GGCGGA", "GGCGTA", "GGCTAA", "GGCTCA", "GGCTGA", "GGCTTA", "GTCAAA",
               "GTCACA", "GTCAGA", "GTCATA", "GTCCAA", "GTCCCA", "GTCCGA", "GTCCTA",
               "GTCGAA", "GTCGCA", "GTCGGA", "GTCGTA", "GTCTAA", "GTCTCA", "GTCTGA",
               "GTCTTA", "TACAAA", "TACACA", "TACAGA", "TACATA", "TACCAA", "TACCCA",
               "TACCGA", "TACCTA", "TACGAA", "TACGCA", "TACGGA", "TACGTA", "TACTAA",
               "TACTCA", "TACTGA", "TACTTA", "TCCAAA", "TCCACA", "TCCAGA", "TCCATA",
               "TCCCAA", "TCCCCA", "TCCCGA", "TCCCTA", "TCCGAA", "TCCGCA", "TCCGGA",
               "TCCGTA", "TCCTAA", "TCCTCA", "TCCTGA", "TCCTTA", "TGCAAA", "TGCACA",
               "TGCAGA", "TGCATA", "TGCCAA", "TGCCCA", "TGCCGA", "TGCCTA", "TGCGAA",
               "TGCGCA", "TGCGGA", "TGCGTA", "TGCTAA", "TGCTCA", "TGCTGA", "TGCTTA",
               "TTCAAA", "TTCACA", "TTCAGA", "TTCATA", "TTCCAA", "TTCCCA", "TTCCGA",
               "TTCCTA", "TTCGAA", "TTCGCA", "TTCGGA", "TTCGTA", "TTCTAA", "TTCTCA",
               "TTCTGA", "TTCTTA", "AACAAG", "AACACG", "AACAGG", "AACATG", "AACCAG",
               "AACCCG", "AACCGG", "AACCTG", "AACGAG", "AACGCG", "AACGGG", "AACGTG",
               "AACTAG", "AACTCG", "AACTGG", "AACTTG", "ACCAAG", "ACCACG", "ACCAGG",
               "ACCATG", "ACCCAG", "ACCCCG", "ACCCGG", "ACCCTG", "ACCGAG", "ACCGCG",
               "ACCGGG", "ACCGTG", "ACCTAG", "ACCTCG", "ACCTGG", "ACCTTG", "AGCAAG",
               "AGCACG", "AGCAGG", "AGCATG", "AGCCAG", "AGCCCG", "AGCCGG", "AGCCTG",
               "AGCGAG", "AGCGCG", "AGCGGG", "AGCGTG", "AGCTAG", "AGCTCG", "AGCTGG",
               "AGCTTG", "ATCAAG", "ATCACG", "ATCAGG", "ATCATG", "ATCCAG", "ATCCCG",
               "ATCCGG", "ATCCTG", "ATCGAG", "ATCGCG", "ATCGGG", "ATCGTG", "ATCTAG",
               "ATCTCG", "ATCTGG", "ATCTTG", "CACAAG", "CACACG", "CACAGG", "CACATG",
               "CACCAG", "CACCCG", "CACCGG", "CACCTG", "CACGAG", "CACGCG", "CACGGG",
               "CACGTG", "CACTAG", "CACTCG", "CACTGG", "CACTTG", "CCCAAG", "CCCACG",
               "CCCAGG", "CCCATG", "CCCCAG", "CCCCCG", "CCCCGG", "CCCCTG", "CCCGAG",
               "CCCGCG", "CCCGGG", "CCCGTG", "CCCTAG", "CCCTCG", "CCCTGG", "CCCTTG",
               "CGCAAG", "CGCACG", "CGCAGG", "CGCATG", "CGCCAG", "CGCCCG", "CGCCGG",
               "CGCCTG", "CGCGAG", "CGCGCG", "CGCGGG", "CGCGTG", "CGCTAG", "CGCTCG",
               "CGCTGG", "CGCTTG", "CTCAAG", "CTCACG", "CTCAGG", "CTCATG", "CTCCAG",
               "CTCCCG", "CTCCGG", "CTCCTG", "CTCGAG", "CTCGCG", "CTCGGG", "CTCGTG",
               "CTCTAG", "CTCTCG", "CTCTGG", "CTCTTG", "GACAAG", "GACACG", "GACAGG",
               "GACATG", "GACCAG", "GACCCG", "GACCGG", "GACCTG", "GACGAG", "GACGCG",
               "GACGGG", "GACGTG", "GACTAG", "GACTCG", "GACTGG", "GACTTG", "GCCAAG",
               "GCCACG", "GCCAGG", "GCCATG", "GCCCAG", "GCCCCG", "GCCCGG", "GCCCTG",
               "GCCGAG", "GCCGCG", "GCCGGG", "GCCGTG", "GCCTAG", "GCCTCG", "GCCTGG",
               "GCCTTG", "GGCAAG", "GGCACG", "GGCAGG", "GGCATG", "GGCCAG", "GGCCCG",
               "GGCCGG", "GGCCTG", "GGCGAG", "GGCGCG", "GGCGGG", "GGCGTG", "GGCTAG",
               "GGCTCG", "GGCTGG", "GGCTTG", "GTCAAG", "GTCACG", "GTCAGG", "GTCATG",
               "GTCCAG", "GTCCCG", "GTCCGG", "GTCCTG", "GTCGAG", "GTCGCG", "GTCGGG",
               "GTCGTG", "GTCTAG", "GTCTCG", "GTCTGG", "GTCTTG", "TACAAG", "TACACG",
               "TACAGG", "TACATG", "TACCAG", "TACCCG", "TACCGG", "TACCTG", "TACGAG",
               "TACGCG", "TACGGG", "TACGTG", "TACTAG", "TACTCG", "TACTGG", "TACTTG",
               "TCCAAG", "TCCACG", "TCCAGG", "TCCATG", "TCCCAG", "TCCCCG", "TCCCGG",
               "TCCCTG", "TCCGAG", "TCCGCG", "TCCGGG", "TCCGTG", "TCCTAG", "TCCTCG",
               "TCCTGG", "TCCTTG", "TGCAAG", "TGCACG", "TGCAGG", "TGCATG", "TGCCAG",
               "TGCCCG", "TGCCGG", "TGCCTG", "TGCGAG", "TGCGCG", "TGCGGG", "TGCGTG",
               "TGCTAG", "TGCTCG", "TGCTGG", "TGCTTG", "TTCAAG", "TTCACG", "TTCAGG",
               "TTCATG", "TTCCAG", "TTCCCG", "TTCCGG", "TTCCTG", "TTCGAG", "TTCGCG",
               "TTCGGG", "TTCGTG", "TTCTAG", "TTCTCG", "TTCTGG", "TTCTTG", "AACAAT",
               "AACACT", "AACAGT", "AACATT", "AACCAT", "AACCCT", "AACCGT", "AACCTT",
               "AACGAT", "AACGCT", "AACGGT", "AACGTT", "AACTAT", "AACTCT", "AACTGT",
               "AACTTT", "ACCAAT", "ACCACT", "ACCAGT", "ACCATT", "ACCCAT", "ACCCCT",
               "ACCCGT", "ACCCTT", "ACCGAT", "ACCGCT", "ACCGGT", "ACCGTT", "ACCTAT",
               "ACCTCT", "ACCTGT", "ACCTTT", "AGCAAT", "AGCACT", "AGCAGT", "AGCATT",
               "AGCCAT", "AGCCCT", "AGCCGT", "AGCCTT", "AGCGAT", "AGCGCT", "AGCGGT",
               "AGCGTT", "AGCTAT", "AGCTCT", "AGCTGT", "AGCTTT", "ATCAAT", "ATCACT",
               "ATCAGT", "ATCATT", "ATCCAT", "ATCCCT", "ATCCGT", "ATCCTT", "ATCGAT",
               "ATCGCT", "ATCGGT", "ATCGTT", "ATCTAT", "ATCTCT", "ATCTGT", "ATCTTT",
               "CACAAT", "CACACT", "CACAGT", "CACATT", "CACCAT", "CACCCT", "CACCGT",
               "CACCTT", "CACGAT", "CACGCT", "CACGGT", "CACGTT", "CACTAT", "CACTCT",
               "CACTGT", "CACTTT", "CCCAAT", "CCCACT", "CCCAGT", "CCCATT", "CCCCAT",
               "CCCCCT", "CCCCGT", "CCCCTT", "CCCGAT", "CCCGCT", "CCCGGT", "CCCGTT",
               "CCCTAT", "CCCTCT", "CCCTGT", "CCCTTT", "CGCAAT", "CGCACT", "CGCAGT",
               "CGCATT", "CGCCAT", "CGCCCT", "CGCCGT", "CGCCTT", "CGCGAT", "CGCGCT",
               "CGCGGT", "CGCGTT", "CGCTAT", "CGCTCT", "CGCTGT", "CGCTTT", "CTCAAT",
               "CTCACT", "CTCAGT", "CTCATT", "CTCCAT", "CTCCCT", "CTCCGT", "CTCCTT",
               "CTCGAT", "CTCGCT", "CTCGGT", "CTCGTT", "CTCTAT", "CTCTCT", "CTCTGT",
               "CTCTTT", "GACAAT", "GACACT", "GACAGT", "GACATT", "GACCAT", "GACCCT",
               "GACCGT", "GACCTT", "GACGAT", "GACGCT", "GACGGT", "GACGTT", "GACTAT",
               "GACTCT", "GACTGT", "GACTTT", "GCCAAT", "GCCACT", "GCCAGT", "GCCATT",
               "GCCCAT", "GCCCCT", "GCCCGT", "GCCCTT", "GCCGAT", "GCCGCT", "GCCGGT",
               "GCCGTT", "GCCTAT", "GCCTCT", "GCCTGT", "GCCTTT", "GGCAAT", "GGCACT",
               "GGCAGT", "GGCATT", "GGCCAT", "GGCCCT", "GGCCGT", "GGCCTT", "GGCGAT",
               "GGCGCT", "GGCGGT", "GGCGTT", "GGCTAT", "GGCTCT", "GGCTGT", "GGCTTT",
               "GTCAAT", "GTCACT", "GTCAGT", "GTCATT", "GTCCAT", "GTCCCT", "GTCCGT",
               "GTCCTT", "GTCGAT", "GTCGCT", "GTCGGT", "GTCGTT", "GTCTAT", "GTCTCT",
               "GTCTGT", "GTCTTT", "TACAAT", "TACACT", "TACAGT", "TACATT", "TACCAT",
               "TACCCT", "TACCGT", "TACCTT", "TACGAT", "TACGCT", "TACGGT", "TACGTT",
               "TACTAT", "TACTCT", "TACTGT", "TACTTT", "TCCAAT", "TCCACT", "TCCAGT",
               "TCCATT", "TCCCAT", "TCCCCT", "TCCCGT", "TCCCTT", "TCCGAT", "TCCGCT",
               "TCCGGT", "TCCGTT", "TCCTAT", "TCCTCT", "TCCTGT", "TCCTTT", "TGCAAT",
               "TGCACT", "TGCAGT", "TGCATT", "TGCCAT", "TGCCCT", "TGCCGT", "TGCCTT",
               "TGCGAT", "TGCGCT", "TGCGGT", "TGCGTT", "TGCTAT", "TGCTCT", "TGCTGT",
               "TGCTTT", "TTCAAT", "TTCACT", "TTCAGT", "TTCATT", "TTCCAT", "TTCCCT",
               "TTCCGT", "TTCCTT", "TTCGAT", "TTCGCT", "TTCGGT", "TTCGTT", "TTCTAT",
               "TTCTCT", "TTCTGT", "TTCTTT", "AATAAA", "AATACA", "AATAGA", "AATATA",
               "AATCAA", "AATCCA", "AATCGA", "AATCTA", "AATGAA", "AATGCA", "AATGGA",
               "AATGTA", "AATTAA", "AATTCA", "AATTGA", "AATTTA", "ACTAAA", "ACTACA",
               "ACTAGA", "ACTATA", "ACTCAA", "ACTCCA", "ACTCGA", "ACTCTA", "ACTGAA",
               "ACTGCA", "ACTGGA", "ACTGTA", "ACTTAA", "ACTTCA", "ACTTGA", "ACTTTA",
               "AGTAAA", "AGTACA", "AGTAGA", "AGTATA", "AGTCAA", "AGTCCA", "AGTCGA",
               "AGTCTA", "AGTGAA", "AGTGCA", "AGTGGA", "AGTGTA", "AGTTAA", "AGTTCA",
               "AGTTGA", "AGTTTA", "ATTAAA", "ATTACA", "ATTAGA", "ATTATA", "ATTCAA",
               "ATTCCA", "ATTCGA", "ATTCTA", "ATTGAA", "ATTGCA", "ATTGGA", "ATTGTA",
               "ATTTAA", "ATTTCA", "ATTTGA", "ATTTTA", "CATAAA", "CATACA", "CATAGA",
               "CATATA", "CATCAA", "CATCCA", "CATCGA", "CATCTA", "CATGAA", "CATGCA",
               "CATGGA", "CATGTA", "CATTAA", "CATTCA", "CATTGA", "CATTTA", "CCTAAA",
               "CCTACA", "CCTAGA", "CCTATA", "CCTCAA", "CCTCCA", "CCTCGA", "CCTCTA",
               "CCTGAA", "CCTGCA", "CCTGGA", "CCTGTA", "CCTTAA", "CCTTCA", "CCTTGA",
               "CCTTTA", "CGTAAA", "CGTACA", "CGTAGA", "CGTATA", "CGTCAA", "CGTCCA",
               "CGTCGA", "CGTCTA", "CGTGAA", "CGTGCA", "CGTGGA", "CGTGTA", "CGTTAA",
               "CGTTCA", "CGTTGA", "CGTTTA", "CTTAAA", "CTTACA", "CTTAGA", "CTTATA",
               "CTTCAA", "CTTCCA", "CTTCGA", "CTTCTA", "CTTGAA", "CTTGCA", "CTTGGA",
               "CTTGTA", "CTTTAA", "CTTTCA", "CTTTGA", "CTTTTA", "GATAAA", "GATACA",
               "GATAGA", "GATATA", "GATCAA", "GATCCA", "GATCGA", "GATCTA", "GATGAA",
               "GATGCA", "GATGGA", "GATGTA", "GATTAA", "GATTCA", "GATTGA", "GATTTA",
               "GCTAAA", "GCTACA", "GCTAGA", "GCTATA", "GCTCAA", "GCTCCA", "GCTCGA",
               "GCTCTA", "GCTGAA", "GCTGCA", "GCTGGA", "GCTGTA", "GCTTAA", "GCTTCA",
               "GCTTGA", "GCTTTA", "GGTAAA", "GGTACA", "GGTAGA", "GGTATA", "GGTCAA",
               "GGTCCA", "GGTCGA", "GGTCTA", "GGTGAA", "GGTGCA", "GGTGGA", "GGTGTA",
               "GGTTAA", "GGTTCA", "GGTTGA", "GGTTTA", "GTTAAA", "GTTACA", "GTTAGA",
               "GTTATA", "GTTCAA", "GTTCCA", "GTTCGA", "GTTCTA", "GTTGAA", "GTTGCA",
               "GTTGGA", "GTTGTA", "GTTTAA", "GTTTCA", "GTTTGA", "GTTTTA", "TATAAA",
               "TATACA", "TATAGA", "TATATA", "TATCAA", "TATCCA", "TATCGA", "TATCTA",
               "TATGAA", "TATGCA", "TATGGA", "TATGTA", "TATTAA", "TATTCA", "TATTGA",
               "TATTTA", "TCTAAA", "TCTACA", "TCTAGA", "TCTATA", "TCTCAA", "TCTCCA",
               "TCTCGA", "TCTCTA", "TCTGAA", "TCTGCA", "TCTGGA", "TCTGTA", "TCTTAA",
               "TCTTCA", "TCTTGA", "TCTTTA", "TGTAAA", "TGTACA", "TGTAGA", "TGTATA",
               "TGTCAA", "TGTCCA", "TGTCGA", "TGTCTA", "TGTGAA", "TGTGCA", "TGTGGA",
               "TGTGTA", "TGTTAA", "TGTTCA", "TGTTGA", "TGTTTA", "TTTAAA", "TTTACA",
               "TTTAGA", "TTTATA", "TTTCAA", "TTTCCA", "TTTCGA", "TTTCTA", "TTTGAA",
               "TTTGCA", "TTTGGA", "TTTGTA", "TTTTAA", "TTTTCA", "TTTTGA", "TTTTTA",
               "AATAAC", "AATACC", "AATAGC", "AATATC", "AATCAC", "AATCCC", "AATCGC",
               "AATCTC", "AATGAC", "AATGCC", "AATGGC", "AATGTC", "AATTAC", "AATTCC",
               "AATTGC", "AATTTC", "ACTAAC", "ACTACC", "ACTAGC", "ACTATC", "ACTCAC",
               "ACTCCC", "ACTCGC", "ACTCTC", "ACTGAC", "ACTGCC", "ACTGGC", "ACTGTC",
               "ACTTAC", "ACTTCC", "ACTTGC", "ACTTTC", "AGTAAC", "AGTACC", "AGTAGC",
               "AGTATC", "AGTCAC", "AGTCCC", "AGTCGC", "AGTCTC", "AGTGAC", "AGTGCC",
               "AGTGGC", "AGTGTC", "AGTTAC", "AGTTCC", "AGTTGC", "AGTTTC", "ATTAAC",
               "ATTACC", "ATTAGC", "ATTATC", "ATTCAC", "ATTCCC", "ATTCGC", "ATTCTC",
               "ATTGAC", "ATTGCC", "ATTGGC", "ATTGTC", "ATTTAC", "ATTTCC", "ATTTGC",
               "ATTTTC", "CATAAC", "CATACC", "CATAGC", "CATATC", "CATCAC", "CATCCC",
               "CATCGC", "CATCTC", "CATGAC", "CATGCC", "CATGGC", "CATGTC", "CATTAC",
               "CATTCC", "CATTGC", "CATTTC", "CCTAAC", "CCTACC", "CCTAGC", "CCTATC",
               "CCTCAC", "CCTCCC", "CCTCGC", "CCTCTC", "CCTGAC", "CCTGCC", "CCTGGC",
               "CCTGTC", "CCTTAC", "CCTTCC", "CCTTGC", "CCTTTC", "CGTAAC", "CGTACC",
               "CGTAGC", "CGTATC", "CGTCAC", "CGTCCC", "CGTCGC", "CGTCTC", "CGTGAC",
               "CGTGCC", "CGTGGC", "CGTGTC", "CGTTAC", "CGTTCC", "CGTTGC", "CGTTTC",
               "CTTAAC", "CTTACC", "CTTAGC", "CTTATC", "CTTCAC", "CTTCCC", "CTTCGC",
               "CTTCTC", "CTTGAC", "CTTGCC", "CTTGGC", "CTTGTC", "CTTTAC", "CTTTCC",
               "CTTTGC", "CTTTTC", "GATAAC", "GATACC", "GATAGC", "GATATC", "GATCAC",
               "GATCCC", "GATCGC", "GATCTC", "GATGAC", "GATGCC", "GATGGC", "GATGTC",
               "GATTAC", "GATTCC", "GATTGC", "GATTTC", "GCTAAC", "GCTACC", "GCTAGC",
               "GCTATC", "GCTCAC", "GCTCCC", "GCTCGC", "GCTCTC", "GCTGAC", "GCTGCC",
               "GCTGGC", "GCTGTC", "GCTTAC", "GCTTCC", "GCTTGC", "GCTTTC", "GGTAAC",
               "GGTACC", "GGTAGC", "GGTATC", "GGTCAC", "GGTCCC", "GGTCGC", "GGTCTC",
               "GGTGAC", "GGTGCC", "GGTGGC", "GGTGTC", "GGTTAC", "GGTTCC", "GGTTGC",
               "GGTTTC", "GTTAAC", "GTTACC", "GTTAGC", "GTTATC", "GTTCAC", "GTTCCC",
               "GTTCGC", "GTTCTC", "GTTGAC", "GTTGCC", "GTTGGC", "GTTGTC", "GTTTAC",
               "GTTTCC", "GTTTGC", "GTTTTC", "TATAAC", "TATACC", "TATAGC", "TATATC",
               "TATCAC", "TATCCC", "TATCGC", "TATCTC", "TATGAC", "TATGCC", "TATGGC",
               "TATGTC", "TATTAC", "TATTCC", "TATTGC", "TATTTC", "TCTAAC", "TCTACC",
               "TCTAGC", "TCTATC", "TCTCAC", "TCTCCC", "TCTCGC", "TCTCTC", "TCTGAC",
               "TCTGCC", "TCTGGC", "TCTGTC", "TCTTAC", "TCTTCC", "TCTTGC", "TCTTTC",
               "TGTAAC", "TGTACC", "TGTAGC", "TGTATC", "TGTCAC", "TGTCCC", "TGTCGC",
               "TGTCTC", "TGTGAC", "TGTGCC", "TGTGGC", "TGTGTC", "TGTTAC", "TGTTCC",
               "TGTTGC", "TGTTTC", "TTTAAC", "TTTACC", "TTTAGC", "TTTATC", "TTTCAC",
               "TTTCCC", "TTTCGC", "TTTCTC", "TTTGAC", "TTTGCC", "TTTGGC", "TTTGTC",
               "TTTTAC", "TTTTCC", "TTTTGC", "TTTTTC", "AATAAG", "AATACG", "AATAGG",
               "AATATG", "AATCAG", "AATCCG", "AATCGG", "AATCTG", "AATGAG", "AATGCG",
               "AATGGG", "AATGTG", "AATTAG", "AATTCG", "AATTGG", "AATTTG", "ACTAAG",
               "ACTACG", "ACTAGG", "ACTATG", "ACTCAG", "ACTCCG", "ACTCGG", "ACTCTG",
               "ACTGAG", "ACTGCG", "ACTGGG", "ACTGTG", "ACTTAG", "ACTTCG", "ACTTGG",
               "ACTTTG", "AGTAAG", "AGTACG", "AGTAGG", "AGTATG", "AGTCAG", "AGTCCG",
               "AGTCGG", "AGTCTG", "AGTGAG", "AGTGCG", "AGTGGG", "AGTGTG", "AGTTAG",
               "AGTTCG", "AGTTGG", "AGTTTG", "ATTAAG", "ATTACG", "ATTAGG", "ATTATG",
               "ATTCAG", "ATTCCG", "ATTCGG", "ATTCTG", "ATTGAG", "ATTGCG", "ATTGGG",
               "ATTGTG", "ATTTAG", "ATTTCG", "ATTTGG", "ATTTTG", "CATAAG", "CATACG",
               "CATAGG", "CATATG", "CATCAG", "CATCCG", "CATCGG", "CATCTG", "CATGAG",
               "CATGCG", "CATGGG", "CATGTG", "CATTAG", "CATTCG", "CATTGG", "CATTTG",
               "CCTAAG", "CCTACG", "CCTAGG", "CCTATG", "CCTCAG", "CCTCCG", "CCTCGG",
               "CCTCTG", "CCTGAG", "CCTGCG", "CCTGGG", "CCTGTG", "CCTTAG", "CCTTCG",
               "CCTTGG", "CCTTTG", "CGTAAG", "CGTACG", "CGTAGG", "CGTATG", "CGTCAG",
               "CGTCCG", "CGTCGG", "CGTCTG", "CGTGAG", "CGTGCG", "CGTGGG", "CGTGTG",
               "CGTTAG", "CGTTCG", "CGTTGG", "CGTTTG", "CTTAAG", "CTTACG", "CTTAGG",
               "CTTATG", "CTTCAG", "CTTCCG", "CTTCGG", "CTTCTG", "CTTGAG", "CTTGCG",
               "CTTGGG", "CTTGTG", "CTTTAG", "CTTTCG", "CTTTGG", "CTTTTG", "GATAAG",
               "GATACG", "GATAGG", "GATATG", "GATCAG", "GATCCG", "GATCGG", "GATCTG",
               "GATGAG", "GATGCG", "GATGGG", "GATGTG", "GATTAG", "GATTCG", "GATTGG",
               "GATTTG", "GCTAAG", "GCTACG", "GCTAGG", "GCTATG", "GCTCAG", "GCTCCG",
               "GCTCGG", "GCTCTG", "GCTGAG", "GCTGCG", "GCTGGG", "GCTGTG", "GCTTAG",
               "GCTTCG", "GCTTGG", "GCTTTG", "GGTAAG", "GGTACG", "GGTAGG", "GGTATG",
               "GGTCAG", "GGTCCG", "GGTCGG", "GGTCTG", "GGTGAG", "GGTGCG", "GGTGGG",
               "GGTGTG", "GGTTAG", "GGTTCG", "GGTTGG", "GGTTTG", "GTTAAG", "GTTACG",
               "GTTAGG", "GTTATG", "GTTCAG", "GTTCCG", "GTTCGG", "GTTCTG", "GTTGAG",
               "GTTGCG", "GTTGGG", "GTTGTG", "GTTTAG", "GTTTCG", "GTTTGG", "GTTTTG",
               "TATAAG", "TATACG", "TATAGG", "TATATG", "TATCAG", "TATCCG", "TATCGG",
               "TATCTG", "TATGAG", "TATGCG", "TATGGG", "TATGTG", "TATTAG", "TATTCG",
               "TATTGG", "TATTTG", "TCTAAG", "TCTACG", "TCTAGG", "TCTATG", "TCTCAG",
               "TCTCCG", "TCTCGG", "TCTCTG", "TCTGAG", "TCTGCG", "TCTGGG", "TCTGTG",
               "TCTTAG", "TCTTCG", "TCTTGG", "TCTTTG", "TGTAAG", "TGTACG", "TGTAGG",
               "TGTATG", "TGTCAG", "TGTCCG", "TGTCGG", "TGTCTG", "TGTGAG", "TGTGCG",
               "TGTGGG", "TGTGTG", "TGTTAG", "TGTTCG", "TGTTGG", "TGTTTG", "TTTAAG",
               "TTTACG", "TTTAGG", "TTTATG", "TTTCAG", "TTTCCG", "TTTCGG", "TTTCTG",
               "TTTGAG", "TTTGCG", "TTTGGG", "TTTGTG", "TTTTAG", "TTTTCG", "TTTTGG",
               "TTTTTG"), NULL))

  vcf.df <- ReadStrelkaVCF("data-raw/cis_8wks_05_cl4_SNVresult.vcf")
  stopifnot(nrow(vcf.df) ==  3544)
  SNS.vcf <- SplitStrelkaSNSVCF(vcf.df)$SNS.vcf

  SNS.vcf <- AddSequence(SNS.vcf)
  CheckSeqContextInVCF(SNS.vcf, "seq.21context")

  SNS.vcf <- AddTranscript(SNS.vcf, .trans.ranges)

  cats <- CreateOneColSNSCatalog(SNS.vcf)

  stopifnot(sum(cats$cat96) == 3336)
  stopifnot(sum(cats$cat1536) == 3336)
  stopifnot(sum(cats$cat192) == 1520)

  stopifnot(cats$cat96 == expected.cat96)

  stopifnot(cats$cat192 == expected.cat192)

  stopifnot(cats$cat1536 == expected.cat1536)

  cat("ok\n")
}

#' This function is to test whether the predefined functions
#' are working correctly to produce the desired DNS catalogs from Strelka VCF.
#' @export
TestStrelkaDNSCatalog <- function() {
  expected.cat.78 <-
    structure(
      c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
        0L, 0L, 10L, 3L, 3L, 0L, 0L, 0L, 2L, 1L, 1L, 2L, 0L, 0L, 0L,
        0L, 0L, 20L, 11L, 4L, 0L, 0L, 0L, 3L, 0L, 0L, 2L, 0L, 0L, 1L,
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 2L, 7L, 0L, 0L, 0L, 0L, 0L,
        0L, 0L, 1L, 5L, 1L, 0L, 2L, 0L, 0L, 15L, 0L, 0L, 1L, 0L, 0L,
        0L, 0L, 0L, 0L),
      .Dim = c(78L, 1L),
      .Dimnames = list(
        c("ACCA",
          "ACCG", "ACCT", "ACGA", "ACGG", "ACGT", "ACTA", "ACTG", "ACTT",
          "ATCA", "ATCC", "ATCG", "ATGA", "ATGC", "ATTA", "CCAA", "CCAG",
          "CCAT", "CCGA", "CCGG", "CCGT", "CCTA", "CCTG", "CCTT", "CGAT",
          "CGGC", "CGGT", "CGTA", "CGTC", "CGTT", "CTAA", "CTAC", "CTAG",
          "CTGA", "CTGC", "CTGG", "CTTA", "CTTC", "CTTG", "GCAA", "GCAG",
          "GCAT", "GCCA", "GCCG", "GCTA", "TAAT", "TACG", "TACT", "TAGC",
          "TAGG", "TAGT", "TCAA", "TCAG", "TCAT", "TCCA", "TCCG", "TCCT",
          "TCGA", "TCGG", "TCGT", "TGAA", "TGAC", "TGAT", "TGCA", "TGCC",
          "TGCT", "TGGA", "TGGC", "TGGT", "TTAA", "TTAC", "TTAG", "TTCA",
          "TTCC", "TTCG", "TTGA", "TTGC", "TTGG"), NULL))

  expected.cat.144 <-
    structure(
      c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 4L, 1L, 0L, 2L, 0L, 0L,
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 1L, 0L, 0L, 0L,
        2L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 2L, 0L, 0L, 0L, 0L,
        0L, 0L, 6L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
        1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 1L,
        1L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L,
        0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 4L, 0L, 0L, 0L, 0L, 0L, 0L,
        0L, 0L, 0L),
      .Dim = c(144L, 1L),
      .Dimnames = list(
        c("AACC",
          "AACG", "AACT", "AAGC", "AAGG", "AAGT", "AATC", "AATG", "AATT",
          "ACCA", "ACCG", "ACCT", "ACGA", "ACGG", "ACGT", "ACTA", "ACTG",
          "ACTT", "AGCA", "AGCC", "AGCT", "AGGA", "AGGC", "AGGT", "AGTA",
          "AGTC", "AGTT", "ATCA", "ATCC", "ATCG", "ATGA", "ATGC", "ATGG",
          "ATTA", "ATTC", "ATTG", "CAAC", "CAAG", "CAAT", "CAGC", "CAGG",
          "CAGT", "CATC", "CATG", "CATT", "CCAA", "CCAG", "CCAT", "CCGA",
          "CCGG", "CCGT", "CCTA", "CCTG", "CCTT", "CGAA", "CGAC", "CGAT",
          "CGGA", "CGGC", "CGGT", "CGTA", "CGTC", "CGTT", "CTAA", "CTAC",
          "CTAG", "CTGA", "CTGC", "CTGG", "CTTA", "CTTC", "CTTG", "GAAC",
          "GAAG", "GAAT", "GACC", "GACG", "GACT", "GATC", "GATG", "GATT",
          "GCAA", "GCAG", "GCAT", "GCCA", "GCCG", "GCCT", "GCTA", "GCTG",
          "GCTT", "GGAA", "GGAC", "GGAT", "GGCA", "GGCC", "GGCT", "GGTA",
          "GGTC", "GGTT", "GTAA", "GTAC", "GTAG", "GTCA", "GTCC", "GTCG",
          "GTTA", "GTTC", "GTTG", "TAAC", "TAAG", "TAAT", "TACC", "TACG",
          "TACT", "TAGC", "TAGG", "TAGT", "TCAA", "TCAG", "TCAT", "TCCA",
          "TCCG", "TCCT", "TCGA", "TCGG", "TCGT", "TGAA", "TGAC", "TGAT",
          "TGCA", "TGCC", "TGCT", "TGGA", "TGGC", "TGGT", "TTAA", "TTAC",
          "TTAG", "TTCA", "TTCC", "TTCG", "TTGA", "TTGC", "TTGG"), NULL))

  expected.cat.136 <-
    structure(
      c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L,
        0L, 1L, 0L, 7L, 1L, 3L, 6L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 2L,
        0L, 0L, 2L, 1L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
        0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 2L, 2L, 5L, 4L, 0L, 0L, 0L,
        0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 5L, 1L, 0L, 0L, 0L, 0L, 0L,
        0L, 0L, 0L, 0L, 0L, 3L, 3L, 0L, 8L, 0L, 0L, 2L, 2L, 1L, 0L, 0L,
        1L, 0L, 0L, 1L, 0L, 0L, 3L, 0L, 0L, 8L, 0L, 0L, 0L, 0L, 0L, 0L,
        0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L,
        0L, 0L, 3L, 0L, 0L, 4L, 0L, 0L, 0L, 1L, 0L),
      .Dim = c(136L, 1L),
      .Dimnames = list(
        c("AACA", "AACC", "AACG", "AACT", "AATA",
          "AATC", "AATG", "AATT", "ACCA", "ACCC", "ACCG", "ACCT", "ACGA",
          "ACGC", "ACGG", "ACGT", "ACTA", "ACTC", "ACTG", "ACTT", "AGCA",
          "AGCC", "AGCG", "AGCT", "ATAA", "ATAC", "ATAG", "ATAT", "ATCA",
          "ATCC", "ATCG", "ATCT", "ATGA", "ATGC", "ATGG", "ATGT", "ATTA",
          "ATTC", "ATTG", "ATTT", "CACA", "CACC", "CACG", "CACT", "CATA",
          "CATC", "CATG", "CCCA", "CCCC", "CCCG", "CCCT", "CCGA", "CCGC",
          "CCGG", "CCTA", "CCTC", "CCTG", "CCTT", "CGCA", "CGCC", "CGCG",
          "CTAA", "CTAC", "CTAG", "CTCA", "CTCC", "CTCG", "CTCT", "CTGA",
          "CTGC", "CTGG", "CTGT", "CTTA", "CTTC", "CTTG", "CTTT", "GACA",
          "GACC", "GACG", "GACT", "GATA", "GATC", "GCCA", "GCCC", "GCCG",
          "GCCT", "GCGA", "GCGC", "GCTA", "GCTC", "GCTG", "GCTT", "GGCA",
          "GGCC", "GTAA", "GTAC", "GTCA", "GTCC", "GTCG", "GTCT", "GTGA",
          "GTGC", "GTGG", "GTGT", "GTTA", "GTTC", "GTTG", "GTTT", "TACA",
          "TACC", "TACG", "TACT", "TATA", "TCCA", "TCCC", "TCCG", "TCCT",
          "TCGA", "TCTA", "TCTC", "TCTG", "TCTT", "TGCA", "TTAA", "TTCA",
          "TTCC", "TTCG", "TTCT", "TTGA", "TTGC", "TTGG", "TTGT", "TTTA",
          "TTTC", "TTTG", "TTTT"), NULL))
  vcf.df <- ReadStrelkaVCF("data-raw/cis_8wks_05_cl4_SNVresult.vcf")
  stopifnot(nrow(vcf.df) ==  3544)
  DNS.vcf <- SplitStrelkaSNSVCF(vcf.df)$DNS.vcf

  DNS.vcf <- AddSequence(DNS.vcf)
  CheckSeqContextInVCF(DNS.vcf, "seq.21context")

  DNS.vcf <- AddTranscript(DNS.vcf, .trans.ranges)

  DNS.cat <- CreateOneColDNSCatalog(DNS.vcf)

  stopifnot(DNS.cat$catDNS78 == expected.cat.78)
  stopifnot(DNS.cat$catDNS144 == expected.cat.144)
  stopifnot(DNS.cat$catQUAD136 == expected.cat.136)

  cat("ok\n")
}

#' This function is to test whether the predefined functions
#' are working correctly to produce the desired SNS and DNS catalogs.
#' @export
TestSNSandDNSCat <- function() {
  # TODO(steve): move this into TestSNSCatalog and TESTSNSCatalog
  # so that these functions are self contained.
  vcf.df <- ReadStrelkaVCF("data-raw/cis_8wks_05_cl4_SNVresult.vcf")
  stopifnot(nrow(vcf.df) ==  3544)
  three.vcfs.df <- SplitSNSVCF(vcf.df)

  # Use default transcript ranges
  TestStrelkaSNSCatalog(three.vcfs.df$SNS.vcf)
  TestStrelakDNSCatalog(three.vcfs.df$DNS.vcf)
}

#' This function is to make catalogs from the sample VCF files
#' to compare with the expected catalog information.
#' @export
TestMakeCatalogFromStrelkaSNSVCFs <- function() {
  # This function is to make catalogs from the sample VCF files
  # to compare with the expected catalog information.

  files <- c("data-raw/HepG2_Cis_1_SNVresult_rmDup.vcf",
             "data-raw/HepG2_Cis_2_SNVresult_rmDup.vcf",
             "data-raw/HepG2_Cis_3_SNVresult_rmDup.vcf",
             "data-raw/HepG2_Cis_4_SNVresult_rmDup.vcf")

  cats <-
    StrelkaVCFFilesToCatalog(
      files,
      genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
      # Use default transcript ranges
      trans.ranges = .trans.ranges)

  prev.catalog.192 <-
    ReadCat192(
      "tests/testthat/testdata/regress.cat.192.csv")
  stopifnot(cats$cat192 == prev.catalog.192)

  prev.catalog.96 <-
    ReadCat96(
      "tests/testthat/testdata/regress.cat.96.csv")
  stopifnot(cats$cat96 == prev.catalog.96)

  prev.catalog.1536 <-
    ReadCat1536(
      "tests/testthat/testdata/regress.cat.1536.csv")
  stopifnot(cats$cat1536 == prev.catalog.1536)

  prev.catalog.DNS.78<-
    ReadCatDNS78(
      "tests/testthat/testdata/regress.cat.dns.78.csv")
  stopifnot(cats$catDNS78 == prev.catalog.DNS.78)

  prev.catalog.QUAD.136<-
    ReadCatQUAD136(
      "tests/testthat/testdata/regress.cat.quad.136.csv")
  stopifnot(cats$catQUAD136 == prev.catalog.QUAD.136)

  prev.catalog.DNS.144<-
    ReadCatDNS144(
      "tests/testthat/testdata/regress.cat.dns.144.csv")
  stopifnot(cats$catDNS144 == prev.catalog.DNS.144)

  cat("ok\n")

  invisible(cats)
}

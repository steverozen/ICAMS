#' This function is to test whether the predefined functions
#' are working correctly to produce the desired SNS catalogs from Strelka VCF.
#'
#' @keywords internal
TestStrelkaVCFToSNSCatalog <- function() {
  expected.catSNS96 <-
    structure(
      c(98L, 247L, 7L, 116L, 184L, 89L, 11L, 232L, 267L,
        570L, 8L, 127L, 109L, 125L, 8L, 167L, 45L, 31L, 5L, 72L, 30L,
        35L, 7L, 66L, 41L, 310L, 5L, 120L, 97L, 41L, 11L, 221L, 56L,
        33L, 19L, 59L, 205L, 583L, 63L, 611L, 150L, 169L, 23L, 118L,
        105L, 134L, 32L, 199L, 68L, 81L, 46L, 63L, 295L, 224L, 333L,
        320L, 8L, 61L, 21L, 22L, 71L, 30L, 31L, 50L, 116L, 57L, 55L,
        74L, 95L, 67L, 106L, 207L, 43L, 134L, 34L, 38L, 81L, 65L, 23L,
        44L, 43L, 25L, 46L, 44L, 23L, 33L, 63L, 60L, 13L, 51L, 28L, 24L,
        47L, 26L, 60L, 42L), .Dim = c(96L, 1L),
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

  expected.catSNS192 <-
    structure(
      c(6L, 4L, 11L, 7L, 12L, 5L, 17L, 10L, 6L, 17L, 13L,
        6L, 11L, 1L, 4L, 5L, 4L, 8L, 30L, 12L, 3L, 11L, 25L, 9L, 11L,
        23L, 17L, 14L, 17L, 10L, 30L, 28L, 15L, 2L, 66L, 14L, 8L, 4L,
        96L, 8L, 8L, 16L, 70L, 25L, 9L, 1L, 70L, 11L, 10L, 42L, 1L, 19L,
        36L, 14L, 3L, 44L, 47L, 125L, 0L, 26L, 25L, 21L, 2L, 30L, 7L,
        4L, 2L, 15L, 6L, 7L, 4L, 13L, 7L, 71L, 0L, 29L, 24L, 6L, 3L,
        42L, 8L, 7L, 4L, 6L, 35L, 79L, 5L, 74L, 38L, 40L, 7L, 11L, 20L,
        32L, 7L, 44L, 36L, 24L, 177L, 16L, 8L, 4L, 15L, 5L, 25L, 24L,
        162L, 9L, 21L, 28L, 47L, 12L, 42L, 18L, 15L, 18L, 2L, 3L, 1L,
        2L, 15L, 66L, 9L, 14L, 19L, 9L, 8L, 9L, 31L, 32L, 44L, 22L, 1L,
        1L, 2L, 1L, 26L, 126L, 20L, 57L, 20L, 49L, 33L, 24L, 16L, 8L,
        7L, 14L, 46L, 30L, 56L, 43L, 1L, 9L, 2L, 4L, 15L, 5L, 2L, 5L,
        16L, 8L, 4L, 18L, 15L, 11L, 16L, 55L, 9L, 30L, 9L, 6L, 10L, 11L,
        4L, 12L, 11L, 4L, 10L, 9L, 2L, 4L, 13L, 15L, 2L, 8L, 7L, 6L,
        5L, 5L, 10L, 6L),
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

  expected.catSNS1536 <-
    structure(
      c(14L, 6L, 5L, 8L, 16L, 37L, 3L, 53L, 0L, 1L, 0L, 0L,
        9L, 10L, 13L, 8L, 13L, 8L, 7L, 8L, 6L, 8L, 4L, 16L, 0L, 0L, 0L,
        0L, 9L, 13L, 21L, 27L, 43L, 29L, 36L, 46L, 91L, 59L, 7L, 122L,
        0L, 0L, 0L, 0L, 9L, 7L, 11L, 22L, 8L, 5L, 9L, 10L, 12L, 4L, 0L,
        13L, 1L, 1L, 0L, 1L, 12L, 8L, 14L, 10L, 9L, 3L, 6L, 7L, 9L, 26L,
        0L, 30L, 0L, 0L, 0L, 0L, 5L, 9L, 9L, 11L, 16L, 12L, 24L, 10L,
        13L, 0L, 0L, 4L, 0L, 1L, 3L, 2L, 13L, 10L, 14L, 25L, 0L, 2L,
        1L, 1L, 8L, 8L, 1L, 27L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 2L, 7L,
        5L, 6L, 6L, 3L, 8L, 1L, 31L, 0L, 2L, 0L, 0L, 7L, 1L, 13L, 16L,
        8L, 3L, 3L, 4L, 3L, 8L, 0L, 9L, 0L, 1L, 2L, 0L, 2L, 2L, 1L, 9L,
        5L, 1L, 6L, 5L, 6L, 1L, 2L, 4L, 1L, 0L, 0L, 0L, 4L, 7L, 11L,
        9L, 24L, 7L, 14L, 9L, 12L, 17L, 0L, 29L, 0L, 1L, 1L, 1L, 11L,
        4L, 13L, 9L, 4L, 4L, 3L, 5L, 6L, 3L, 0L, 10L, 1L, 0L, 0L, 0L,
        10L, 8L, 4L, 9L, 7L, 7L, 3L, 5L, 6L, 23L, 1L, 23L, 0L, 1L, 1L,
        1L, 7L, 7L, 9L, 5L, 16L, 14L, 15L, 24L, 11L, 7L, 1L, 6L, 1L,
        1L, 0L, 2L, 11L, 7L, 14L, 37L, 16L, 9L, 20L, 10L, 48L, 54L, 2L,
        85L, 2L, 0L, 2L, 1L, 10L, 3L, 7L, 16L, 10L, 8L, 7L, 12L, 10L,
        11L, 1L, 12L, 0L, 0L, 1L, 1L, 12L, 12L, 8L, 23L, 4L, 4L, 4L,
        5L, 7L, 1L, 0L, 4L, 0L, 0L, 0L, 1L, 6L, 4L, 3L, 5L, 0L, 1L, 1L,
        0L, 1L, 3L, 0L, 4L, 1L, 0L, 1L, 1L, 4L, 2L, 6L, 3L, 8L, 4L, 7L,
        5L, 88L, 45L, 5L, 75L, 0L, 0L, 1L, 1L, 12L, 20L, 21L, 23L, 10L,
        3L, 6L, 4L, 2L, 1L, 0L, 4L, 1L, 0L, 0L, 0L, 11L, 7L, 36L, 19L,
        6L, 1L, 6L, 1L, 3L, 4L, 0L, 6L, 1L, 0L, 1L, 0L, 5L, 5L, 4L, 11L,
        4L, 2L, 5L, 3L, 4L, 0L, 0L, 2L, 0L, 0L, 0L, 1L, 3L, 7L, 4L, 10L,
        0L, 0L, 1L, 0L, 2L, 1L, 0L, 4L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 2L,
        15L, 4L, 16L, 4L, 5L, 3L, 0L, 4L, 5L, 0L, 1L, 0L, 12L, 8L, 26L,
        15L, 2L, 1L, 4L, 4L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 5L,
        5L, 5L, 4L, 0L, 3L, 2L, 2L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 2L, 3L,
        4L, 1L, 1L, 2L, 1L, 1L, 27L, 11L, 2L, 17L, 0L, 0L, 2L, 0L, 3L,
        3L, 8L, 8L, 4L, 0L, 5L, 2L, 5L, 3L, 0L, 1L, 1L, 0L, 1L, 1L, 6L,
        2L, 15L, 3L, 0L, 0L, 3L, 0L, 2L, 2L, 1L, 0L, 0L, 0L, 0L, 2L,
        2L, 3L, 5L, 3L, 1L, 1L, 1L, 2L, 10L, 3L, 0L, 4L, 0L, 0L, 1L,
        1L, 0L, 5L, 6L, 6L, 2L, 1L, 6L, 2L, 17L, 7L, 1L, 8L, 0L, 0L,
        0L, 1L, 7L, 3L, 3L, 5L, 13L, 1L, 4L, 6L, 3L, 5L, 0L, 5L, 1L,
        0L, 0L, 0L, 11L, 6L, 25L, 19L, 4L, 3L, 3L, 11L, 2L, 5L, 0L, 3L,
        1L, 1L, 0L, 2L, 2L, 5L, 5L, 5L, 4L, 4L, 9L, 7L, 20L, 10L, 2L,
        14L, 0L, 2L, 1L, 2L, 7L, 19L, 20L, 22L, 20L, 12L, 17L, 18L, 31L,
        15L, 4L, 41L, 2L, 1L, 2L, 1L, 7L, 8L, 4L, 18L, 4L, 1L, 7L, 6L,
        2L, 7L, 1L, 10L, 2L, 1L, 2L, 3L, 5L, 8L, 8L, 16L, 2L, 3L, 3L,
        4L, 2L, 2L, 1L, 3L, 1L, 5L, 1L, 2L, 5L, 3L, 7L, 4L, 10L, 13L,
        20L, 18L, 69L, 0L, 13L, 61L, 3L, 6L, 3L, 5L, 12L, 47L, 16L, 59L,
        0L, 1L, 3L, 1L, 1L, 3L, 2L, 2L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 3L,
        3L, 8L, 14L, 9L, 9L, 12L, 2L, 26L, 5L, 1L, 3L, 4L, 7L, 8L, 17L,
        27L, 2L, 1L, 2L, 4L, 3L, 2L, 0L, 2L, 1L, 2L, 0L, 0L, 2L, 2L,
        3L, 3L, 4L, 2L, 6L, 3L, 5L, 2L, 0L, 2L, 2L, 3L, 1L, 2L, 2L, 8L,
        3L, 4L, 5L, 4L, 11L, 7L, 6L, 2L, 1L, 11L, 1L, 4L, 3L, 2L, 4L,
        3L, 6L, 11L, 0L, 3L, 5L, 3L, 3L, 6L, 0L, 14L, 0L, 0L, 1L, 2L,
        3L, 3L, 2L, 22L, 3L, 1L, 3L, 7L, 5L, 3L, 0L, 0L, 0L, 0L, 1L,
        2L, 5L, 3L, 1L, 4L, 16L, 36L, 20L, 33L, 140L, 98L, 13L, 134L,
        6L, 6L, 7L, 14L, 72L, 108L, 96L, 116L, 10L, 10L, 13L, 18L, 14L,
        10L, 1L, 25L, 3L, 2L, 1L, 0L, 12L, 9L, 12L, 20L, 4L, 6L, 20L,
        12L, 10L, 16L, 3L, 13L, 2L, 0L, 3L, 3L, 9L, 18L, 19L, 27L, 10L,
        5L, 10L, 4L, 4L, 10L, 1L, 34L, 2L, 2L, 2L, 7L, 11L, 5L, 6L, 4L,
        49L, 30L, 58L, 53L, 39L, 22L, 5L, 28L, 55L, 35L, 47L, 39L, 49L,
        13L, 47L, 39L, 2L, 1L, 1L, 0L, 3L, 3L, 0L, 12L, 2L, 2L, 3L, 2L,
        0L, 6L, 1L, 2L, 4L, 1L, 2L, 4L, 3L, 1L, 0L, 4L, 1L, 1L, 1L, 2L,
        3L, 2L, 1L, 8L, 7L, 4L, 3L, 1L, 1L, 3L, 0L, 7L, 2L, 4L, 2L, 2L,
        6L, 1L, 3L, 2L, 12L, 8L, 13L, 9L, 19L, 26L, 4L, 17L, 13L, 9L,
        14L, 18L, 13L, 11L, 19L, 16L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 3L,
        1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 6L, 4L, 4L, 2L, 1L, 0L, 4L,
        4L, 2L, 1L, 4L, 3L, 3L, 7L, 4L, 1L, 1L, 5L, 4L, 2L, 2L, 0L, 6L,
        1L, 0L, 4L, 0L, 3L, 3L, 4L, 1L, 7L, 2L, 2L, 7L, 4L, 5L, 0L, 4L,
        9L, 9L, 13L, 8L, 7L, 7L, 7L, 8L, 0L, 0L, 0L, 1L, 3L, 1L, 1L,
        5L, 1L, 1L, 2L, 1L, 0L, 1L, 0L, 2L, 5L, 4L, 0L, 2L, 1L, 1L, 2L,
        2L, 0L, 1L, 1L, 1L, 2L, 4L, 1L, 2L, 3L, 2L, 4L, 4L, 1L, 3L, 0L,
        7L, 6L, 5L, 1L, 6L, 4L, 2L, 2L, 6L, 14L, 9L, 13L, 9L, 13L, 15L,
        3L, 20L, 11L, 14L, 15L, 24L, 10L, 15L, 25L, 34L, 2L, 1L, 0L,
        0L, 8L, 0L, 0L, 20L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 4L, 8L, 4L,
        3L, 16L, 2L, 2L, 0L, 5L, 4L, 3L, 2L, 3L, 1L, 4L, 5L, 0L, 10L,
        9L, 10L, 8L, 4L, 5L, 1L, 15L, 5L, 8L, 4L, 4L, 3L, 5L, 2L, 4L,
        9L, 3L, 3L, 10L, 1L, 2L, 0L, 3L, 6L, 5L, 10L, 7L, 29L, 10L, 15L,
        2L, 2L, 4L, 5L, 4L, 24L, 4L, 1L, 12L, 2L, 2L, 3L, 5L, 1L, 0L,
        1L, 4L, 5L, 3L, 7L, 6L, 0L, 1L, 0L, 5L, 1L, 0L, 0L, 5L, 2L, 3L,
        2L, 4L, 4L, 7L, 5L, 6L, 2L, 3L, 0L, 8L, 2L, 3L, 1L, 4L, 0L, 5L,
        6L, 11L, 9L, 2L, 2L, 4L, 7L, 9L, 0L, 5L, 5L, 5L, 5L, 7L, 6L,
        7L, 11L, 8L, 0L, 0L, 0L, 1L, 3L, 1L, 0L, 3L, 0L, 1L, 0L, 0L,
        0L, 1L, 0L, 0L, 4L, 1L, 2L, 2L, 1L, 1L, 0L, 6L, 2L, 0L, 0L, 2L,
        3L, 2L, 4L, 2L, 6L, 0L, 3L, 5L, 0L, 4L, 0L, 5L, 0L, 0L, 5L, 1L,
        3L, 4L, 4L, 3L, 0L, 0L, 4L, 4L, 1L, 3L, 0L, 1L, 4L, 1L, 4L, 2L,
        18L, 5L, 12L, 4L, 0L, 1L, 2L, 5L, 16L, 1L, 1L, 9L, 1L, 2L, 2L,
        1L, 1L, 3L, 4L, 1L, 1L, 1L, 4L, 0L, 4L, 3L, 0L, 6L, 1L, 0L, 0L,
        0L, 1L, 0L, 1L, 3L, 15L, 8L, 12L, 8L, 3L, 3L, 0L, 4L, 4L, 4L,
        7L, 3L, 3L, 4L, 6L, 11L, 8L, 5L, 17L, 15L, 10L, 16L, 0L, 9L,
        7L, 12L, 16L, 10L, 24L, 16L, 12L, 28L, 5L, 2L, 6L, 6L, 31L, 6L,
        10L, 12L, 4L, 1L, 6L, 4L, 7L, 4L, 8L, 3L, 11L, 8L, 10L, 16L,
        9L, 10L, 1L, 18L, 4L, 1L, 2L, 5L, 6L, 5L, 6L, 0L, 7L, 3L, 4L,
        3L, 3L, 0L, 0L, 1L, 2L, 2L, 5L, 6L, 4L, 2L, 3L, 7L, 2L, 2L, 2L,
        3L, 2L, 0L, 1L, 2L, 6L, 4L, 5L, 5L, 6L, 1L, 5L, 12L, 2L, 0L,
        2L, 0L, 8L, 0L, 0L, 13L, 3L, 0L, 1L, 3L, 3L, 4L, 1L, 4L, 5L,
        2L, 0L, 5L, 6L, 1L, 0L, 6L, 4L, 0L, 1L, 0L, 4L, 4L, 3L, 5L, 4L,
        1L, 2L, 1L, 4L, 2L, 0L, 5L, 1L, 6L, 5L, 0L, 6L, 2L, 3L, 4L, 3L,
        0L, 2L, 1L, 4L, 4L, 0L, 1L, 1L, 6L, 6L, 1L, 1L, 2L, 1L, 7L, 1L,
        0L, 1L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 6L,
        0L, 2L, 4L, 3L, 0L, 0L, 2L, 6L, 4L, 3L, 3L, 5L, 1L, 4L, 2L, 5L,
        0L, 2L, 1L, 0L, 0L, 0L, 1L, 3L, 1L, 6L, 0L, 2L, 0L, 2L, 2L, 1L,
        2L, 0L, 1L, 4L, 3L, 0L, 0L, 0L, 4L, 7L, 4L, 3L, 2L, 2L, 3L, 2L,
        0L, 0L, 0L, 1L, 1L, 1L, 8L, 1L, 0L, 6L, 2L, 1L, 2L, 1L, 1L, 2L,
        1L, 1L, 1L, 3L, 0L, 0L, 1L, 6L, 1L, 1L, 1L, 3L, 1L, 2L, 2L, 4L,
        0L, 0L, 6L, 1L, 1L, 1L, 6L, 1L, 0L, 6L, 2L, 1L, 0L, 1L, 5L, 2L,
        0L, 0L, 2L, 5L, 2L, 0L, 5L, 5L, 4L, 2L, 3L, 4L, 1L, 3L, 7L, 3L,
        0L, 1L, 1L, 4L, 1L, 1L, 11L, 0L, 1L, 8L, 2L, 2L, 1L, 1L, 2L,
        8L, 0L, 3L, 7L, 1L, 1L, 0L, 2L, 7L, 10L, 4L, 9L, 3L, 1L, 2L,
        0L),
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

  vcf.df <-
    ReadStrelkaSNSVCF(system.file("extdata",
                               "MCF10A_Carb_Low_cl2_Strelka_SNS.vcf",
                               package = "ICAMS",
                               mustWork = TRUE))
  stopifnot(nrow(vcf.df) ==  10293)
  SNS.vcf <- SplitStrelkaSNSVCF(vcf.df)$SNS.vcf

  SNS.vcf <- AddSequence(SNS.vcf, genome = "GRCh37")
  CheckSeqContextInVCF(SNS.vcf, "seq.21context")

  SNS.vcf <- AddTranscript(SNS.vcf, trans.ranges.GRCh37)

  cats <- CreateOneColSNSCatalog(SNS.vcf)

  stopifnot(sum(cats$catSNS96) == 9652)
  stopifnot(sum(cats$catSNS1536) == 9652)
  stopifnot(sum(cats$catSNS192) == 3878)

  stopifnot(cats$catSNS96 == expected.catSNS96)

  stopifnot(cats$catSNS192 == expected.catSNS192)

  stopifnot(cats$catSNS1536 == expected.catSNS1536)

  cat("ok\n")
}

#' This function is to test whether the predefined functions
#' are working correctly to produce the desired DNS catalogs from Strelka VCF.
#' @keywords internal
TestStrelkaVCFToDNSCatalog <- function() {
  expected.cat.78 <-
    structure(
      c(0L, 0L, 0L, 0L, 0L, 0L, 5L, 0L, 2L, 0L, 0L, 0L, 0L,
        0L, 0L, 22L, 7L, 12L, 0L, 1L, 1L, 2L, 1L, 4L, 0L, 0L, 1L, 0L,
        0L, 2L, 55L, 33L, 3L, 4L, 2L, 0L, 3L, 2L, 0L, 14L, 5L, 5L, 1L,
        0L, 2L, 2L, 0L, 2L, 0L, 0L, 1L, 28L, 10L, 30L, 2L, 0L, 2L, 1L,
        0L, 2L, 0L, 2L, 11L, 1L, 2L, 4L, 1L, 0L, 19L, 2L, 1L, 0L, 1L,
        0L, 0L, 0L, 0L, 0L),
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
        0L, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 1L, 0L, 2L, 0L, 1L, 5L, 0L, 0L,
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 2L, 4L, 0L, 0L, 1L, 0L, 0L, 0L,
        3L, 0L, 3L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L,
        0L, 0L, 12L, 4L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 8L, 0L,
        0L, 4L, 1L, 1L, 10L, 2L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 3L, 2L,
        1L, 2L, 0L, 0L, 4L, 1L, 0L, 4L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L,
        0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 1L, 2L, 0L, 0L, 0L,
        0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 7L, 1L, 1L, 0L, 1L,
        0L, 0L, 0L, 0L, 0L),
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
      c(1L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 2L, 0L, 3L, 0L,
        0L, 1L, 0L, 8L, 3L, 18L, 13L, 5L, 12L, 0L, 7L, 0L, 0L, 1L, 1L,
        7L, 9L, 0L, 13L, 1L, 0L, 15L, 0L, 0L, 0L, 1L, 0L, 2L, 0L, 0L,
        0L, 0L, 0L, 0L, 2L, 0L, 2L, 2L, 1L, 0L, 0L, 8L, 12L, 5L, 14L,
        0L, 0L, 0L, 2L, 1L, 0L, 0L, 5L, 1L, 7L, 2L, 4L, 6L, 1L, 3L, 0L,
        0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 9L, 0L, 24L, 0L, 0L, 4L,
        0L, 4L, 2L, 2L, 0L, 0L, 0L, 4L, 4L, 0L, 19L, 0L, 0L, 2L, 1L,
        0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 4L, 1L, 2L, 1L,
        7L, 1L, 1L, 0L, 0L, 0L, 0L, 6L, 2L, 1L, 5L, 0L, 0L, 0L, 0L, 0L),
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
  vcf.df <-
    ReadStrelkaSNSVCF(system.file("extdata",
                               "MCF10A_Carb_Low_cl2_Strelka_SNS.vcf",
                               package = "ICAMS",
                               mustWork = TRUE))
  stopifnot(nrow(vcf.df) ==  10293)
  DNS.vcf <- SplitStrelkaSNSVCF(vcf.df)$DNS.vcf

  DNS.vcf <- AddSequence(DNS.vcf, genome = "GRCh37")
  CheckSeqContextInVCF(DNS.vcf, "seq.21context")

  DNS.vcf <- AddTranscript(DNS.vcf, trans.ranges.GRCh37)

  DNS.cat <- CreateOneColDNSCatalog(DNS.vcf)

  stopifnot(sum(DNS.cat$catDNS78) == 313)
  stopifnot(sum(DNS.cat$catDNS144) == 113)
  stopifnot(sum(DNS.cat$catDNS136) == 313)

  stopifnot(DNS.cat$catDNS78 == expected.cat.78)
  stopifnot(DNS.cat$catDNS144 == expected.cat.144)
  stopifnot(DNS.cat$catDNS136 == expected.cat.136)

  cat("ok\n")
}

#' This function is to make catalogs from the sample Mutect VCF file
#' to compare with the expected catalog information.
#' @keywords internal
TestMakeCatalogFromMutectVCFs <- function() {
  files <- c(system.file("extdata",
                         "MCF10A_Carb_Low_cl2_Mutect.vcf",
                         package = "ICAMS",
                         mustWork = TRUE))

  cats <-
    MutectVCFFilesToCatalog(
      files,
      genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
      # Use default transcript ranges
      trans.ranges = trans.ranges.GRCh37)

  prev.catalog.192 <-
    ReadCatSNS192(
      system.file("extdata",
                  "mutect.regress.cat.sns.192.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catSNS192 == prev.catalog.192)

  prev.catalog.96 <-
    ReadCatSNS96(
      system.file("extdata",
                  "mutect.regress.cat.sns.96.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catSNS96 == prev.catalog.96)

  prev.catalog.1536 <-
    ReadCatSNS1536(
      system.file("extdata",
                  "mutect.regress.cat.sns.1536.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catSNS1536 == prev.catalog.1536)

  prev.catalog.DNS.78<-
    ReadCatDNS78(
      system.file("extdata",
                  "mutect.regress.cat.dns.78.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catDNS78 == prev.catalog.DNS.78)

  prev.catalog.DNS.136<-
    ReadCatDNS136(
      system.file("extdata",
                  "mutect.regress.cat.dns.136.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catQUAD136 == prev.catalog.DNS.136)

  prev.catalog.DNS.144<-
    ReadCatDNS144(
      system.file("extdata",
                  "mutect.regress.cat.dns.144.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catDNS144 == prev.catalog.DNS.144)

  prev.catalog.indels<-
    ReadCatID(
      system.file("extdata",
                  "mutect.regress.cat.indels.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catID == prev.catalog.indels)

  cat("ok\n")

  invisible(cats)
}

#' This function is to make catalogs from the sample Strelka SNS VCF files
#' to compare with the expected catalog information.
#' @keywords internal
TestMakeCatalogFromStrelkaSNSVCFs <- function() {
  files <- c(system.file("extdata",
                         "MCF10A_Carb_Low_cl2_Strelka_SNS.vcf",
                         package = "ICAMS",
                         mustWork = TRUE))

  cats <-
    StrelkaSNSVCFFilesToCatalog(
      files,
      genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
      # Use default transcript ranges
      trans.ranges = trans.ranges.GRCh37)

  prev.catalog.192 <-
    ReadCatSNS192(
      system.file("extdata",
                  "new.regress.cat.sns.192.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catSNS192 == prev.catalog.192)

  prev.catalog.96 <-
    ReadCatSNS96(
      system.file("extdata",
                  "new.regress.cat.sns.96.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catSNS96 == prev.catalog.96)

  prev.catalog.1536 <-
    ReadCatSNS1536(
      system.file("extdata",
                  "new.regress.cat.sns.1536.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catSNS1536 == prev.catalog.1536)

  prev.catalog.DNS.78<-
    ReadCatDNS78(
      system.file("extdata",
                  "new.regress.cat.dns.78.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catDNS78 == prev.catalog.DNS.78)

  prev.catalog.DNS.136<-
    ReadCatDNS136(
      system.file("extdata",
                  "new.regress.cat.dns.136.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catQUAD136 == prev.catalog.DNS.136)

  prev.catalog.DNS.144<-
    ReadCatDNS144(
      system.file("extdata",
                  "new.regress.cat.dns.144.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cats$catDNS144 == prev.catalog.DNS.144)

  cat("ok\n")

  invisible(cats)
}

#' This function is to make catalogs from the sample Strelka ID VCF files
#' to compare with the expected catalog information.
#' @keywords internal
TestMakeCatalogFromStrelkaIDVCFs <- function() {
  files <- c(system.file("extdata",
                         "MCF10A_Carb_Low_cl2_INDELresult.vcf",
                         package = "ICAMS",
                         mustWork = TRUE))

  cat.ID <-
    StrelkaIDVCFFilesToCatalog(
      files,
      genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5)

  prev.catalog.indels<-
    ReadCatID(
      system.file("extdata",
                  "new.regress.cat.indels.csv",
                  package = "ICAMS",
                  mustWork = TRUE))
  stopifnot(cat.ID == prev.catalog.indels)

  cat("ok\n")

  invisible(cat.ID)
}

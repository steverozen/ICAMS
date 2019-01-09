context("CanonicalizeDNS")

test_that("CanonicalizeDNS function is working properly", {
  to.test <-
    c("AACC", "AACG", "AACT", "AAGC", "AAGG", "AAGT", "AATC", "AATG",
      "AATT", "ACCA", "ACCG", "ACCT", "ACGA", "ACGG", "ACGT", "ACTA",
      "ACTG", "ACTT", "AGCA", "AGCC", "AGCT", "AGGA", "AGGC", "AGGT",
      "AGTA", "AGTC", "AGTT", "ATCA", "ATCC", "ATCG", "ATGA", "ATGC",
      "ATGG", "ATTA", "ATTC", "ATTG", "CAAC", "CAAG", "CAAT", "CAGC",
      "CAGG", "CAGT", "CATC", "CATG", "CATT", "CCAA", "CCAG", "CCAT",
      "CCGA", "CCGG", "CCGT", "CCTA", "CCTG", "CCTT", "CGAA", "CGAC",
      "CGAT", "CGGA", "CGGC", "CGGT", "CGTA", "CGTC", "CGTT", "CTAA",
      "CTAC", "CTAG", "CTGA", "CTGC", "CTGG", "CTTA", "CTTC", "CTTG",
      "GAAC", "GAAG", "GAAT", "GACC", "GACG", "GACT", "GATC", "GATG",
      "GATT", "GCAA", "GCAG", "GCAT", "GCCA", "GCCG", "GCCT", "GCTA",
      "GCTG", "GCTT", "GGAA", "GGAC", "GGAT", "GGCA", "GGCC", "GGCT",
      "GGTA", "GGTC", "GGTT", "GTAA", "GTAC", "GTAG", "GTCA", "GTCC",
      "GTCG", "GTTA", "GTTC", "GTTG", "TAAC", "TAAG", "TAAT", "TACC",
      "TACG", "TACT", "TAGC", "TAGG", "TAGT", "TCAA", "TCAG", "TCAT",
      "TCCA", "TCCG", "TCCT", "TCGA", "TCGG", "TCGT", "TGAA", "TGAC",
      "TGAT", "TGCA", "TGCC", "TGCT", "TGGA", "TGGC", "TGGT", "TTAA",
      "TTAC", "TTAG", "TTCA", "TTCC", "TTCG", "TTGA", "TTGC", "TTGG"
    )
  result <-
    c(AACC = "TTGG", AACG = "TTCG", AACT = "TTAG", AAGC = "TTGC",
      AAGG = "TTCC", AAGT = "TTAC", AATC = "TTGA", AATG = "TTCA", AATT = "TTAA",
      ACCA = "ACCA", ACCG = "ACCG", ACCT = "ACCT", ACGA = "ACGA", ACGG = "ACGG",
      ACGT = "ACGT", ACTA = "ACTA", ACTG = "ACTG", ACTT = "ACTT", AGCA = "CTTG",
      AGCC = "CTGG", AGCT = "CTAG", AGGA = "CTTC", AGGC = "CTGC", AGGT = "CTAC",
      AGTA = "CTTA", AGTC = "CTGA", AGTT = "CTAA", ATCA = "ATCA", ATCC = "ATCC",
      ATCG = "ATCG", ATGA = "ATGA", ATGC = "ATGC", ATGG = "ATCC", ATTA = "ATTA",
      ATTC = "ATGA", ATTG = "ATCA", CAAC = "TGGT", CAAG = "TGCT", CAAT = "TGAT",
      CAGC = "TGGC", CAGG = "TGCC", CAGT = "TGAC", CATC = "TGGA", CATG = "TGCA",
      CATT = "TGAA", CCAA = "CCAA", CCAG = "CCAG", CCAT = "CCAT", CCGA = "CCGA",
      CCGG = "CCGG", CCGT = "CCGT", CCTA = "CCTA", CCTG = "CCTG", CCTT = "CCTT",
      CGAA = "CGTT", CGAC = "CGGT", CGAT = "CGAT", CGGA = "CGTC", CGGC = "CGGC",
      CGGT = "CGGT", CGTA = "CGTA", CGTC = "CGTC", CGTT = "CGTT", CTAA = "CTAA",
      CTAC = "CTAC", CTAG = "CTAG", CTGA = "CTGA", CTGC = "CTGC", CTGG = "CTGG",
      CTTA = "CTTA", CTTC = "CTTC", CTTG = "CTTG", GAAC = "TCGT", GAAG = "TCCT",
      GAAT = "TCAT", GACC = "TCGG", GACG = "TCCG", GACT = "TCAG", GATC = "TCGA",
      GATG = "TCCA", GATT = "TCAA", GCAA = "GCAA", GCAG = "GCAG", GCAT = "GCAT",
      GCCA = "GCCA", GCCG = "GCCG", GCCT = "GCAG", GCTA = "GCTA", GCTG = "GCCA",
      GCTT = "GCAA", GGAA = "CCTT", GGAC = "CCGT", GGAT = "CCAT", GGCA = "CCTG",
      GGCC = "CCGG", GGCT = "CCAG", GGTA = "CCTA", GGTC = "CCGA", GGTT = "CCAA",
      GTAA = "ACTT", GTAC = "ACGT", GTAG = "ACCT", GTCA = "ACTG", GTCC = "ACGG",
      GTCG = "ACCG", GTTA = "ACTA", GTTC = "ACGA", GTTG = "ACCA", TAAC = "TAGT",
      TAAG = "TACT", TAAT = "TAAT", TACC = "TAGG", TACG = "TACG", TACT = "TACT",
      TAGC = "TAGC", TAGG = "TAGG", TAGT = "TAGT", TCAA = "TCAA", TCAG = "TCAG",
      TCAT = "TCAT", TCCA = "TCCA", TCCG = "TCCG", TCCT = "TCCT", TCGA = "TCGA",
      TCGG = "TCGG", TCGT = "TCGT", TGAA = "TGAA", TGAC = "TGAC", TGAT = "TGAT",
      TGCA = "TGCA", TGCC = "TGCC", TGCT = "TGCT", TGGA = "TGGA", TGGC = "TGGC",
      TGGT = "TGGT", TTAA = "TTAA", TTAC = "TTAC", TTAG = "TTAG", TTCA = "TTCA",
      TTCC = "TTCC", TTCG = "TTCG", TTGA = "TTGA", TTGC = "TTGC", TTGG = "TTGG"
    )
  expect_equal(CanonicalizeDNS(substr(to.test, 1, 2),
                               substr(to.test, 3, 4)),
               result)
})

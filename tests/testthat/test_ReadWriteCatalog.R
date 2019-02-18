context("Read and Write Catalog")

test_that("Functions for reading and writing catalog are working properly", {
  read.fn <-
    c(ReadCatSNS96,
      ReadCatSNS192,
      ReadCatSNS1536,
      ReadCatDNS78,
      ReadCatDNS144,
      ReadCatDNS136,
      ReadCatID
    )

  write.fn <-
    c(WriteCatSNS96,
      WriteCatSNS192,
      WriteCatSNS1536,
      WriteCatDNS78,
      WriteCatDNS144,
      WriteCatDNS136,
      WriteCatID)

  fl <-
    c("testdata/BTSG_WGS_PCAWG.sns.96.csv",
      "testdata/BTSG_WGS_PCAWG.sns.192.csv",
      "testdata/BTSG_WGS_PCAWG.sns.1536.csv",
      "testdata/BTSG_WGS_PCAWG.dns.78.csv",
      "testdata/HepG2.dns.144.csv",
      "testdata/HepG2.dns.136.csv",
      "testdata/BTSG_WGS_PCAWG.indels.csv"
    )

  Test1Cat <- function(my.read, my.write, my.file) {
    ct1 <- my.read(my.file)
    my.write(ct1, "tmp.ct.txt")
    ct2 <- my.read("tmp.ct.txt")
    expect_equal(ct1, ct2)
  }

  discard <- mapply(Test1Cat, read.fn, write.fn, fl)
  unlink("tmp.ct.txt")
})

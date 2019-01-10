context("Read and Write Catalog")

test_that("Functions for reading and writing catalog are working properly", {
  read.fn <-
    c(ReadCat96,
      ReadCat192,
      ReadCat1536,
      ReadCatDNS78,
      ReadCatDNS144,
      ReadCatQUAD136,
      ReadCatID
    )

  write.fn <-
    c(WriteCat96,
      WriteCat192,
      WriteCat1536,
      WriteCatDNS78,
      WriteCatDNS144,
      WriteCatQUAD136,
      WriteCatID)

  fl <-
    c("testdata/BTSG_WGS_PCAWG.96.csv",
      "testdata/BTSG_WGS_PCAWG.192.csv",
      "testdata/BTSG_WGS_PCAWG.1536.csv",
      "testdata/BTSG_WGS_PCAWG.dinucs.csv",
      "testdata/HepG2.dinucs.144.csv",
      "testdata/HepG2.quad.136.csv",
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

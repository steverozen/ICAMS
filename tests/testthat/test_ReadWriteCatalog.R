context("Read and Write Catalog")

test_that("Functions for reading and writing catalog are working properly", {
  read.fn <-
    c(ReadCatalog,
      ReadCatalog,
      ReadCatalog,
      ReadCatalog,
      ReadCatalog,
      ReadCatalog,
      ReadCatalog)

  write.fn <-
    c(WriteCatalog,
      WriteCatalog,
      WriteCatalog,
      WriteCatalog,
      WriteCatalog,
      WriteCatalog,
      WriteCatalog)

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
    ct1 <- my.read(my.file, ref.genome = "GRCh37",
                   region = "genome", catalog.type = "counts")
    my.write(ct1, "tmp.ct.txt")
    ct2 <- my.read("tmp.ct.txt", ref.genome = "GRCh37",
                   region = "genome", catalog.type = "counts")
    expect_equal(ct1, ct2)
  }

  discard <- mapply(Test1Cat, read.fn, write.fn, fl)
  unlink("tmp.ct.txt")
})

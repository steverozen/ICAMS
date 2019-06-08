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
    c("testdata/BTSG_WGS_PCAWG.sbs.96.csv",
      "testdata/BTSG_WGS_PCAWG.sbs.192.csv",
      "testdata/BTSG_WGS_PCAWG.sbs.1536.csv",
      "testdata/BTSG_WGS_PCAWG.dbs.78.csv",
      "testdata/HepG2.dbs.144.csv",
      "testdata/HepG2.dbs.136.csv",
      "testdata/BTSG_WGS_PCAWG.indels.csv"
    )

  Test1Cat <- function(my.read, my.write, my.file) {
    ct1 <- my.read(my.file, ref.genome = "GRCh37",
                   region = "genome", catalog.type = "counts")
    my.write(ct1, paste0(tempdir(), "\\tmp.ct.txt"))
    ct2 <- my.read(paste0(tempdir(), "\\tmp.ct.txt"), ref.genome = "GRCh37",
                   region = "genome", catalog.type = "counts")
    expect_equal(ct1, ct2)
  }

  discard <- mapply(Test1Cat, read.fn, write.fn, fl)
  unlink(paste0(tempdir(), "\\tmp.ct.txt"))
})

test_that("Functions for reading and writing catalogs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))

  fl <-
    c("testdata/regress.cat.sbs.96.tsv",
      "testdata/BTSG_WGS_PCAWG.sbs.96.csv",
      "testdata/BTSG_WGS_PCAWG.sbs.192.csv",
      "testdata/BTSG_WGS_PCAWG.sbs.1536.csv",
      "testdata/BTSG_WGS_PCAWG.dbs.78.csv",
      "testdata/HepG2.dbs.144.csv",
      "testdata/HepG2.dbs.136.csv",
      "testdata/BTSG_WGS_PCAWG.indels.csv"
    )
  
  a.region <-
    c("genome",
      "genome",
      "transcript",
      "genome",
      "genome",
      "transcript",
      "genome",
      "genome"
      )
  
  Test1Cat <- function(my.file, my.region) {
    ct1 <- ReadCatalog(my.file, ref.genome = "GRCh37",
                       region = my.region, catalog.type = "counts")
    f <- tempfile("catalog")
    WriteCatalog(ct1, f)
    ct2 <- ReadCatalog(f, ref.genome = "GRCh37",
                       region = my.region, catalog.type = "counts")
    expect_equal(ct1, ct2)
  }
  
  discard <- mapply(Test1Cat, fl, a.region)
  temp.files <- list.files(tempdir(), full.names = TRUE, pattern = "^catalog")
  invisible(file.remove(temp.files)) 
})

test_that("Functions for converting COSMIC and SigProfiler-formatted catalogs to ICAMS", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))

  my.file <-
    c("testdata/COSMIC-Cat/COSMIC_SBS1_SBS96.csv",
      "testdata/COSMIC-Cat/COSMIC_SBS85_SBS96.csv",
      "testdata/COSMIC-Cat/COSMIC_SBS1_TSB192.csv",
      "testdata/COSMIC-Cat/COSMIC_DBS1_DBS78.csv",
      "testdata/COSMIC-Cat/COSMIC_ID1_ID96.csv",
      "testdata/SigPro-Cat/21BRCA.SBS96.tsv",
      "testdata/SigPro-Cat/21BRCA.DBS78.tsv",
      "testdata/SigPro-Cat/Strelka.ID.GRCh37.s1.ID83.tsv",
      "testdata/SigPro-Cat/Strelka.ID.GRCh37.s1.ID96.tsv"
    )
  
  my.region <-
    c("genome",
      "genome",
      "transcript",
      "genome",
      "genome",
      "genome",
      "genome",
      "genome",
      "genome"
      )

  Test1Cat <- function(my.file, my.region) {
    ct1 <- ReadCatalog(file = my.file, ref.genome = "GRCh37",
                       region = my.region, 
                       catalog.type = "counts",
                       stop.on.error = TRUE)
    f <- tempfile("catalogConv")
    WriteCatalog(ct1, f)
    ct2 <- ReadCatalog(f, ref.genome = "GRCh37",
                       region = my.region, catalog.type = "counts")
    expect_equal(ct1, ct2)
  }
  
  discard <- mapply(Test1Cat, my.file, my.region)
  temp.files <- list.files(tempdir(), full.names = TRUE, pattern = "^catalogConv")
  invisible(file.remove(temp.files)) 
})

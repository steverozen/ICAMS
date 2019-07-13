context("AddAndCheckSequenceID")

test_that("AddAndCheckSequenceID function with hg19", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  load("testdata/test_AddAndCheckSequenceID.Rdata")
  id.vcf <- ReadStrelkaIDVCF("testdata/Strelka.ID.GRCh37.vcf")
  df <- 
    AddAndCheckSequenceID(id.vcf, 
                          ref.genome = 
                          BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5)
  df1 <- AddAndCheckSequenceID(id.vcf, ref.genome = "GRCh37")
  df2 <- AddAndCheckSequenceID(id.vcf, ref.genome = "hg19")
  expect_equal(df, strelka.ID.vcf.GRCh37)
  expect_equal(df, df1)
  expect_equal(df, df2)
})

test_that("AddAndCheckSequenceID with hg38", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  load("testdata/test_AddAndCheckSequenceID.Rdata")
  id.vcf <- ReadStrelkaIDVCF("testdata/Strelka.ID.GRCh38.vcf")
  df3 <- AddAndCheckSequenceID(id.vcf, ref.genome = BSgenome.Hsapiens.UCSC.hg38)
  df4 <- AddAndCheckSequenceID(id.vcf, ref.genome = "GRCh38")
  df5 <- AddAndCheckSequenceID(id.vcf, ref.genome = "hg38")
  expect_equal(df3, strelka.ID.vcf.GRCh38)
  expect_equal(df3, df4)
  expect_equal(df3, df5)
})


context("AddAndCheckSequenceID")

test_that("AddAndCheckSequenceID function is working properly", {
  load("testdata/test_AddAndCheckSequenceID.Rdata")
  id.vcf <- ReadStrelkaIDVCF("testdata/Strelka.ID.GRCh37.vcf")
  df <- AddAndCheckSequenceID(id.vcf, genome = BSgenome.Hsapiens.1000genomes.hs37d5)
  df1 <- AddAndCheckSequenceID(id.vcf, genome = "GRCh37")
  df2 <- AddAndCheckSequenceID(id.vcf, genome = "hg19")
  expect_equal(df, strelka.ID.vcf.GRCh37)
  expect_equal(df, df1)
  expect_equal(df, df2)
})

test_that("AddAndCheckSequenceID with hg38", {
  load("testdata/test_AddAndCheckSequenceID.Rdata")
  id.vcf <- ReadStrelkaIDVCF("testdata/Strelka.ID.GRCh38.vcf")
  df3 <- AddAndCheckSequenceID(id.vcf, genome = BSgenome.Hsapiens.UCSC.hg38)
  df4 <- AddAndCheckSequenceID(id.vcf, genome = "GRCh38")
  df5 <- AddAndCheckSequenceID(id.vcf, genome = "hg38")
  expect_equal(df3, strelka.ID.vcf.GRCh38)
  expect_equal(df3, df4)
  expect_equal(df3, df5)

})


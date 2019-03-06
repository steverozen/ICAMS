context("AddAndCheckSequenceID")

test_that("AddAndCheckSequenceID function is working properly", {
  id.vcf <- ReadStrelkaIDVCF("testdata/Strelka.ID.GRCh37.vcf")
  # TODO(steve): review w/ Nanhai, why calling twice?
  strelka.ID.vcf.GRCh37 <-
    AddAndCheckSequenceID(id.vcf, genome = BSgenome.Hsapiens.1000genomes.hs37d5)
  df <- AddAndCheckSequenceID(id.vcf, genome = BSgenome.Hsapiens.1000genomes.hs37d5)
  df1 <- AddAndCheckSequenceID(id.vcf, genome = "GRCh37")
  df2 <- AddAndCheckSequenceID(id.vcf, genome = "hg19")
  expect_equal(df, df1)
  expect_equal(df, df2)
})

test_that("AddAndCheckSequenceID with hg38", {
  id.vcf <- ReadStrelkaIDVCF("testdata/Strelka.ID.GRCh38.vcf")
  strelka.ID.vcf.GRCh38 <-
    AddAndCheckSequenceID(
      id.vcf, genome = BSgenome.Hsapiens.UCSC.hg38,
                          flag.mismatches = 20)

  df3 <- AddAndCheckSequenceID(
    id.vcf, genome = BSgenome.Hsapiens.UCSC.hg38,
    flag.mismatches = 20)

  df4 <- AddAndCheckSequenceID(id.vcf, genome = "GRCh38")
  df5 <- AddAndCheckSequenceID(id.vcf, genome = "hg38")
  expect_equal(df3, df4)
  expect_equal(df3, df5)

})


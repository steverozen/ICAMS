context("AddSequence")

test_that("AddSequence for GRCh37", {
  load("testdata/test_AddSequence.Rdata")
  list.of.vcf <- ReadAndSplitStrelkaSNSVCFs("testdata/Strelka.SNS.GRCh37.vcf")
  sns.vcf <- list.of.vcf$SNS.vcfs[[1]]
  df <- AddSequence(sns.vcf, genome = BSgenome.Hsapiens.1000genomes.hs37d5)
  df1 <- AddSequence(sns.vcf, genome = "GRCh37")
  df2 <- AddSequence(sns.vcf, genome = "hg19")
  expect_equal(df, strelka.SNS.vcf.GRCh37)
  expect_equal(df, df1)
  expect_equal(df, df2)
})

test_that("AddSequence for GRCh38", {

  list.of.vcf <- ReadAndSplitStrelkaSNSVCFs("testdata/Strelka.SNS.GRCh38.vcf")
  sns.vcf <- list.of.vcf$SNS.vcfs[[1]]
  df3 <- AddSequence(sns.vcf, genome = BSgenome.Hsapiens.UCSC.hg38)
  df4 <- AddSequence(sns.vcf, genome = "GRCh38")
  df5 <- AddSequence(sns.vcf, genome = "hg38")
  # expect_equal(df3, strelka.SNS.vcf.GRCh38)
  expect_equal(df3, df4)
  expect_equal(df3, df5)
})

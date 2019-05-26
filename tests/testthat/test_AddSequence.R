context("AddSequence")

test_that("AddSequence for GRCh37", {
  load("testdata/test_AddSequence.Rdata")
  list.of.vcf <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  sns.vcf <- list.of.vcf$SBS.vcfs[[1]]
  df <- AddSequence(sns.vcf, ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5)
  df1 <- AddSequence(sns.vcf, ref.genome = "GRCh37")
  df2 <- AddSequence(sns.vcf, ref.genome = "hg19")
  expect_equal(df, strelka.SBS.vcf.GRCh37)
  expect_equal(df, df1)
  expect_equal(df, df2)
})

test_that("AddSequence for GRCh38", {
  load("testdata/test_AddSequence.Rdata")
  list.of.vcf <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh38.vcf")
  sns.vcf <- list.of.vcf$SBS.vcfs[[1]]
  df3 <- AddSequence(sns.vcf, ref.genome = BSgenome.Hsapiens.UCSC.hg38)
  df4 <- AddSequence(sns.vcf, ref.genome = "GRCh38")
  df5 <- AddSequence(sns.vcf, ref.genome = "hg38")
  expect_equal(df3, strelka.SBS.vcf.GRCh38)
  expect_equal(df3, df4)
  expect_equal(df3, df5)
})

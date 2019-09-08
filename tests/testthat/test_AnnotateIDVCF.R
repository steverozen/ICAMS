context("AnnotateIDVCF")

test_that("AnnotateIDVCF function with hg19", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  load("testdata/test_AnnotateIDVCF.Rdata")
  id.vcf <- ReadStrelkaIDVCF("testdata/Strelka.ID.GRCh37.vcf")
  df <- 
    AnnotateIDVCF(id.vcf, 
                  ref.genome = 
                    BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5)
  df1 <- AnnotateIDVCF(id.vcf, ref.genome = "GRCh37")
  df2 <- AnnotateIDVCF(id.vcf, ref.genome = "hg19")
  expect_equal(df, strelka.ID.vcf.GRCh37)
  expect_equal(df, df1)
  expect_equal(df, df2)
})

test_that("AnnotateIDVCF with hg38", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  load("testdata/test_AnnotateIDVCF.Rdata")
  id.vcf <- ReadStrelkaIDVCF("testdata/Strelka.ID.GRCh38.vcf")
  df3 <- AnnotateIDVCF(id.vcf, 
                       ref.genome = 
                         BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  df4 <- AnnotateIDVCF(id.vcf, ref.genome = "GRCh38")
  df5 <- AnnotateIDVCF(id.vcf, ref.genome = "hg38")
  expect_equal(df3, strelka.ID.vcf.GRCh38)
  expect_equal(df3, df4)
  expect_equal(df3, df5)
})


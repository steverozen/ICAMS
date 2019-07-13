context("AddSeqContext")

test_that("AddSeqContext for GRCh37", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  load("testdata/test_AddSeqContext.Rdata")
  list.of.vcf <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  sbs.vcf <- list.of.vcf$SBS.vcfs[[1]]
  df <- 
    AddSeqContext(sbs.vcf, 
                  ref.genome = 
                  BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5)
  df1 <- AddSeqContext(sbs.vcf, ref.genome = "GRCh37")
  df2 <- AddSeqContext(sbs.vcf, ref.genome = "hg19")
  expect_equal(df, strelka.SBS.vcf.GRCh37)
  expect_equal(df, df1)
  expect_equal(df, df2)
})

test_that("AddSeqContext for GRCh38", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  load("testdata/test_AddSeqContext.Rdata")
  list.of.vcf <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh38.vcf")
  sbs.vcf <- list.of.vcf$SBS.vcfs[[1]]
  df3 <- AddSeqContext(sbs.vcf, ref.genome = BSgenome.Hsapiens.UCSC.hg38)
  df4 <- AddSeqContext(sbs.vcf, ref.genome = "GRCh38")
  df5 <- AddSeqContext(sbs.vcf, ref.genome = "hg38")
  expect_equal(df3, strelka.SBS.vcf.GRCh38)
  expect_equal(df3, df4)
  expect_equal(df3, df5)
})

test_that("AddSeqContext for GRCm38", {
  skip_if("" == system.file(package = "BSgenome.Mmusculus.UCSC.mm10"))
  stopifnot(requireNamespace("BSgenome.Mmusculus.UCSC.mm10"))
  load("testdata/test_AddSeqContext.Rdata")
  list.of.vcf <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCm38.vcf")
  sbs.vcf <- list.of.vcf$SBS.vcfs[[1]]
  df6 <- AddSeqContext(sbs.vcf, ref.genome 
                       = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  df7 <- AddSeqContext(sbs.vcf, ref.genome = "GRCm38")
  df8 <- AddSeqContext(sbs.vcf, ref.genome = "mm10")
  expect_equal(df6, strelka.SBS.vcf.GRCm38)
  expect_equal(df6, df7)
  expect_equal(df6, df8)
})

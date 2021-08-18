context("AnnotateIDVCF")

test_that("AnnotateIDVCF function with hg19", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  load("testdata/test_AnnotateIDVCF.Rdata")
  id.vcf <- ReadStrelkaIDVCF("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf")
  list <- 
    AnnotateIDVCF(id.vcf, 
                  ref.genome = 
                    BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5)
  list1 <- AnnotateIDVCF(id.vcf, ref.genome = "GRCh37")
  list3 <- AnnotateIDVCF(id.vcf, ref.genome = "GRCh37", seq.context.width = 10)
  list2 <- AnnotateIDVCF(id.vcf, ref.genome = "hg19")
  expect_equal(list$annotated.vcf, strelka.ID.vcf.GRCh37)
  expect_equal(list$annotated.vcf, list1$annotated.vcf)
  expect_equal(list$annotated.vcf, list2$annotated.vcf)
})

test_that("AnnotateIDVCF with hg38", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  load("testdata/test_AnnotateIDVCF.Rdata")
  id.vcf <- ReadStrelkaIDVCF("testdata/Strelka.ID.GRCh38.vcf")
  list3 <- AnnotateIDVCF(id.vcf, 
                       ref.genome = 
                         BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  list4 <- AnnotateIDVCF(id.vcf, ref.genome = "GRCh38")
  list5 <- AnnotateIDVCF(id.vcf, ref.genome = "hg38")
  expect_equal(list3$annotated.vcf, strelka.ID.vcf.GRCh38)
  expect_equal(list3$annotated.vcf, list4$annotated.vcf)
  expect_equal(list3$annotated.vcf, list5$annotated.vcf)
})


context("Making catalogs from Strelka GRCh38 VCFs")

test_that("StrelkaSNSVCFFilesToCatalog", {
  cat1 <- StrelkaSNSVCFFilesToCatalog("testdata/Strelka.SNS.GRCh38.vcf",
                                      ref.genome = BSgenome.Hsapiens.UCSC.hg38,
                                      trans.ranges = trans.ranges.GRCh38,
                                      region = "genome")
  cat2 <- StrelkaSNSVCFFilesToCatalog("testdata/Strelka.SNS.GRCh38.vcf",
                                      ref.genome = "GRCh38",
                                      trans.ranges = trans.ranges.GRCh38,
                                      region = "genome")
  cat3 <- StrelkaSNSVCFFilesToCatalog("testdata/Strelka.SNS.GRCh38.vcf",
                                      ref.genome = "hg38",
                                      trans.ranges = trans.ranges.GRCh38,
                                      region = "genome")
  expect_equal(cat1$catSNS96, cat2$catSNS96)
  expect_equal(cat1$catSNS96, cat3$catSNS96)
})

test_that("StrelkaIDVCFFilesToCatalog", {
  cat4 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     ref.genome = BSgenome.Hsapiens.UCSC.hg38,
                                     region = "genome")
  cat5 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     ref.genome = "GRCh38",
                                     region = "genome")
  cat6 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     ref.genome = "hg38",
                                     region = "genome")
  expect_equal(cat4, cat5)
  expect_equal(cat4, cat6)

})

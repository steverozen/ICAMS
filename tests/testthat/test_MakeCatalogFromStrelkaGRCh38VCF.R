context("Making catalogs from Strelka GRCh38 VCFs")

test_that("StrelkaSNSVCFFilesToCatalog", {
  cat1 <- StrelkaSNSVCFFilesToCatalog("testdata/Strelka.SNS.GRCh38.vcf",
                                      genome = BSgenome.Hsapiens.UCSC.hg38,
                                      trans.ranges = trans.ranges.GRCh38)
  cat2 <- StrelkaSNSVCFFilesToCatalog("testdata/Strelka.SNS.GRCh38.vcf",
                                      genome = "GRCh38",
                                      trans.ranges = trans.ranges.GRCh38)
  cat3 <- StrelkaSNSVCFFilesToCatalog("testdata/Strelka.SNS.GRCh38.vcf",
                                      genome = "hg38",
                                      trans.ranges = trans.ranges.GRCh38)
  expect_equal(cat1, cat2)
  expect_equal(cat1, cat3)
})


test_that("StrelkaIDVCFFilesToCatalog", {
  cat4 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                       genome = BSgenome.Hsapiens.UCSC.hg38)
  cat5 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     genome = "GRCh38")
  cat6 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     genome = "hg38")
  expect_equal(cat4, cat5)
  expect_equal(cat4, cat6)

})

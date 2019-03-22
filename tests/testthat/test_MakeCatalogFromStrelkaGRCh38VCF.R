context("Making catalogs from Strelka GRCh38 VCFs")

test_that("StrelkaSNSVCFFilesToCatalog", {
  cat1 <- StrelkaSNSVCFFilesToCatalog("testdata/Strelka.SNS.GRCh38.vcf",
                                      genome = BSgenome.Hsapiens.UCSC.hg38,
                                      trans.ranges = trans.ranges.GRCh38,
                                      region = "genome")
  cat2 <- StrelkaSNSVCFFilesToCatalog("testdata/Strelka.SNS.GRCh38.vcf",
                                      genome = "GRCh38",
                                      trans.ranges = trans.ranges.GRCh38,
                                      region = "genome")
  cat3 <- StrelkaSNSVCFFilesToCatalog("testdata/Strelka.SNS.GRCh38.vcf",
                                      genome = "hg38",
                                      trans.ranges = trans.ranges.GRCh38,
                                      region = "genome")
  expect_equal(cat1$catSNS96$catalog, cat2$catSNS96$catalog)
  expect_equal(cat1$catSNS96$catalog, cat3$catSNS96$catalog)
})


test_that("StrelkaIDVCFFilesToCatalog", {
  cat4 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     genome = BSgenome.Hsapiens.UCSC.hg38,
                                     region = "genome")
  cat5 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     genome = "GRCh38",
                                     region = "genome")
  cat6 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     genome = "hg38",
                                     region = "genome")
  expect_equal(cat4$catSNS96$catalog, cat5$catSNS96$catalog)
  expect_equal(cat4$catSNS96$catalog, cat6$catSNS96$catalog)

})

context("Making catalogs from Strelka GRCh38 VCFs")

test_that("StrelkaSBSVCFFilesToCatalog", {
  cat1 <- StrelkaSBSVCFFilesToCatalog("testdata/Strelka.SBS.GRCh38.vcf",
                                      ref.genome = BSgenome.Hsapiens.UCSC.hg38,
                                      trans.ranges = trans.ranges.GRCh38,
                                      region = "genome")
  cat2 <- StrelkaSBSVCFFilesToCatalog("testdata/Strelka.SBS.GRCh38.vcf",
                                      ref.genome = "GRCh38",
                                      trans.ranges = trans.ranges.GRCh38,
                                      region = "genome")
  cat3 <- StrelkaSBSVCFFilesToCatalog("testdata/Strelka.SBS.GRCh38.vcf",
                                      ref.genome = "hg38",
                                      trans.ranges = trans.ranges.GRCh38,
                                      region = "genome")
  expect_equal(cat1$catSBS96, cat2$catSBS96)
  expect_equal(cat1$catSBS96, cat3$catSBS96)
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

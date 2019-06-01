context("Making catalogs from Mutect GRCh38 VCFs")

test_that("MutectVCFFilesToCatalog", {
  cat1 <- MutectVCFFilesToCatalog("testdata/Mutect.GRCh38.vcf",
                                      ref.genome = BSgenome.Hsapiens.UCSC.hg38,
                                      trans.ranges = trans.ranges.GRCh38,
                                      region = "genome")
  cat2 <- MutectVCFFilesToCatalog("testdata/Mutect.GRCh38.vcf",
                                  ref.genome = "GRCh38",
                                  trans.ranges = trans.ranges.GRCh38,
                                  region = "genome")
  cat3 <- MutectVCFFilesToCatalog("testdata/Mutect.GRCh38.vcf",
                                  ref.genome = "hg38",
                                  trans.ranges = trans.ranges.GRCh38,
                                  region = "genome")
  expect_equal(cat1$catSBS96, cat2$catSBS96)
  expect_equal(cat1$catSBS96, cat3$catSBS96)
})

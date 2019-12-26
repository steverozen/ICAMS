context("Making catalogs from Mutect GRCh38 VCFs")

test_that("MutectVCFFilesToCatalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  cat1 <- MutectVCFFilesToCatalog("testdata/Mutect.GRCh38.vcf",
                                  ref.genome = 
                                  BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
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
  
  # Test for case when trans.ranges is not supplied.
  cat4 <- MutectVCFFilesToCatalog("testdata/Mutect.GRCh38.vcf",
                                  ref.genome = "hg38",
                                  region = "genome")
  cat5 <- MutectVCFFilesToCatalog("testdata/Mutect.GRCh38.vcf",
                                  ref.genome = "hg38")
  expect_equal(cat1$catSBS96, cat2$catSBS96)
  expect_equal(cat1$catSBS96, cat3$catSBS96)
  expect_equal(cat1$catSBS96, cat4$catSBS96)
  expect_equal(attributes(cat5$catSBS96)$region, "unknown")
  expect_null(attributes(cat5$catSBS96)$abundance)
})

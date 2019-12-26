context("Making catalogs from Strelka GRCh38 VCFs")

test_that("StrelkaSBSVCFFilesToCatalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  cat1 <- StrelkaSBSVCFFilesToCatalog("testdata/Strelka.SBS.GRCh38.vcf",
                                      ref.genome = 
                                      BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
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
  cat4 <- StrelkaSBSVCFFilesToCatalog("testdata/Strelka.SBS.GRCh38.vcf",
                                      ref.genome = "hg38",
                                      region = "genome")
  cat5 <- StrelkaSBSVCFFilesToCatalog("testdata/Strelka.SBS.GRCh38.vcf",
                                      ref.genome = "hg38")
  expect_equal(cat1$catSBS96, cat2$catSBS96)
  expect_equal(cat1$catSBS96, cat3$catSBS96)
  expect_equal(cat1$catSBS96, cat4$catSBS96)
  expect_equal(attributes(cat5$catSBS96)$region, "unknown")
  expect_null(attributes(cat5$catSBS96)$abundance)
})

test_that("StrelkaIDVCFFilesToCatalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  cat5 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     ref.genome = 
                                     BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                                     region = "genome")
  cat6 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     ref.genome = "GRCh38",
                                     region = "genome")
  cat7 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     ref.genome = "hg38",
                                     region = "genome")
  cat8 <- StrelkaIDVCFFilesToCatalog("testdata/Strelka.ID.GRCh38.vcf",
                                     ref.genome = "hg38")
  expect_equal(cat5, cat6)
  expect_equal(cat5, cat7)
  expect_equal(attributes(cat8$catalog)$region, "unknown")
  expect_null(attributes(cat8$catalog)$abundance)
})

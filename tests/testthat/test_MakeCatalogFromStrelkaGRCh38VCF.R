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
  cat6 <- VCFsToCatalogs("testdata/Strelka.SBS.GRCh38.vcf", ref.genome = "hg38",
                         variant.caller = "strelka", region = "genome")
  expect_equal(cat1, cat2)
  expect_equal(cat1, cat3)
  expect_equal(cat1, cat4)
  
  cat6$catID <- cat6$catID166 <-NULL
  expect_equal(cat1, cat6)
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
  cat9 <- VCFsToCatalogs("testdata/Strelka.ID.GRCh38.vcf", ref.genome = "hg38",
                         variant.caller = "strelka", region = "genome")
  
  expect_equal(cat5, cat6)
  expect_equal(cat5, cat7)
  expect_equal(cat5$catalog, cat9$catID)
  expect_equal(attributes(cat8$catalog)$region, "unknown")
  expect_null(attributes(cat8$catalog)$abundance)
})

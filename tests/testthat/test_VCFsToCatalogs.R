context("VCFsToCatalogs function")

test_that("VCFsToCatalogs function for Mutect VCFs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Mutect-GRCh37", full.names = TRUE)
  catalogs1 <- MutectVCFFilesToCatalog(files, ref.genome = "hg19", 
                                       region = "genome")
  catalogs2 <- MutectVCFFilesToCatalog(files, ref.genome = "hg19", 
                                       region = "genome",
                                       return.annotated.vcfs = TRUE)
  catalogs3 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "mutect", region = "genome")
  catalogs4 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "mutect", region = "genome",
                              return.annotated.vcfs = TRUE)
  catalogs5 <- VCFsToCatalogs(files, ref.genome = "hg19", region = "genome",
                              get.vaf.function = GetMutectVAF)
  catalogs6 <- VCFsToCatalogs(files, ref.genome = "hg19", region = "genome",
                              get.vaf.function = GetMutectVAF,
                              return.annotated.vcfs = TRUE)
  expect_equal(catalogs1, catalogs3)
  expect_equal(catalogs2, catalogs4)
  expect_equal(catalogs1, catalogs5)
  expect_equal(catalogs2, catalogs6)
})

test_that("VCFsToCatalogs function for Strelka SBS VCFs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Strelka-SBS-GRCh37/", full.names = TRUE)
  catalogs1 <- StrelkaSBSVCFFilesToCatalog(files, ref.genome = "hg19", 
                                           region = "genome")
  catalogs2 <- StrelkaSBSVCFFilesToCatalog(files, ref.genome = "hg19", 
                                           region = "genome",
                                           return.annotated.vcfs = TRUE)
  catalogs3 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "strelka", region = "genome")
  catalogs4 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "strelka", region = "genome",
                              return.annotated.vcfs = TRUE)
  catalogs5 <- VCFsToCatalogs(files, ref.genome = "hg19", region = "genome",
                              get.vaf.function = GetStrelkaVAF)
  catalogs6 <- VCFsToCatalogs(files, ref.genome = "hg19", region = "genome",
                              get.vaf.function = GetStrelkaVAF,
                              return.annotated.vcfs = TRUE)
  catalogs3$catID <- NULL
  expect_equal(catalogs1, catalogs3)
  catalogs5$catID <- NULL
  expect_equal(catalogs1, catalogs5)
  catalogs4$catID <- catalogs4$annotated.vcfs$ID <- NULL
  expect_equal(catalogs2, catalogs4)
  catalogs6$catID <- catalogs6$annotated.vcfs$ID <- NULL
  expect_equal(catalogs2, catalogs6)
})

test_that("VCFsToCatalogs function for Strelka ID VCFs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Strelka-ID-GRCh37/", full.names = TRUE)
  catalogs1 <- StrelkaIDVCFFilesToCatalog(files, ref.genome = "hg19", 
                                           region = "genome")
  catalogs2 <- StrelkaIDVCFFilesToCatalog(files, ref.genome = "hg19", 
                                           region = "genome",
                                           return.annotated.vcfs = TRUE)
  catalogs3 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "strelka", region = "genome")
  catalogs4 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "strelka", region = "genome",
                              return.annotated.vcfs = TRUE)
  expect_equal(catalogs1$catalog, catalogs3$catID)
  expect_equal(catalogs1$discarded.variants, catalogs3$discarded.variants)
  
  expect_equal(catalogs2$catalog, catalogs4$catID)
  expect_equal(catalogs2$discarded.variants, catalogs4$discarded.variants)
  expect_equal(catalogs2$annotated.vcfs, catalogs4$annotated.vcfs$ID)
  unlink("Rplots.pdf")
})
context("VCFsToCatalogs function")

test_that("VCFsToCatalogs function for Mutect VCFs", {
  rlang::with_options(lifecycle_verbosity = "quiet", {
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
                              get.vaf.function = GetMutectVAF,
                              filter.status = "PASS")
  catalogs6 <- VCFsToCatalogs(files, ref.genome = "hg19", region = "genome",
                              get.vaf.function = GetMutectVAF,
                              return.annotated.vcfs = TRUE,
                              filter.status = "PASS")
  expect_equal(catalogs1, catalogs3)
  expect_equal(catalogs2, catalogs4)
  expect_equal(catalogs1, catalogs5)
  expect_equal(catalogs2, catalogs6)
  })
})

test_that("VCFsToCatalogs function for Strelka SBS VCFs", {
  rlang::with_options(lifecycle_verbosity = "quiet", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Strelka-SBS-GRCh37/", full.names = TRUE)
  catalogs1 <- StrelkaSBSVCFFilesToCatalog(files, ref.genome = "hg19", 
                                           region = "genome")
  catalogs2 <- StrelkaSBSVCFFilesToCatalog(files, ref.genome = "hg19", 
                                           region = "genome",
                                           return.annotated.vcfs = TRUE)
  catalogs3 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "strelka", region = "genome",
                              filter.status = NULL)
  catalogs4 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "strelka", region = "genome",
                              return.annotated.vcfs = TRUE,
                              filter.status = NULL)
  catalogs5 <- VCFsToCatalogs(files, ref.genome = "hg19", region = "genome",
                              get.vaf.function = GetStrelkaVAF,
                              filter.status = NULL)
  catalogs6 <- VCFsToCatalogs(files, ref.genome = "hg19", region = "genome",
                              get.vaf.function = GetStrelkaVAF,
                              return.annotated.vcfs = TRUE,
                              filter.status = NULL)
  catalogs3$catID <- catalogs3$catID166 <- NULL
  expect_equal(catalogs1, catalogs3)
  catalogs5$catID <- catalogs5$catID166 <- NULL
  expect_equal(catalogs1, catalogs5)
  catalogs4$catID <- catalogs4$catID166 <- catalogs4$annotated.vcfs$ID <- NULL
  expect_equal(catalogs2, catalogs4)
  catalogs6$catID <- catalogs6$catID166 <- catalogs6$annotated.vcfs$ID <- NULL
  expect_equal(catalogs2, catalogs6)
  })
})

#' @import rlang with_options
test_that("VCFsToCatalogs function for Strelka ID VCFs", {
  rlang::with_options(lifecycle_verbosity = "quiet", {
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
  expect_equal(catalogs1$catID166, catalogs3$catID166)
  expect_equal(catalogs1$discarded.variants, catalogs3$discarded.variants)
  
  expect_equal(catalogs2$catalog, catalogs4$catID)
  expect_equal(catalogs2$catID166, catalogs4$catID166)
  expect_equal(catalogs2$discarded.variants, catalogs4$discarded.variants)
  expect_equal(catalogs2$annotated.vcfs, catalogs4$annotated.vcfs$ID)
  unlink("Rplots.pdf")
  })
})
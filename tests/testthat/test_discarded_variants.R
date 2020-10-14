context("Test functions dealing with discarded variants and annotated vcfs")

test_that("ReadAndSplitMutectVCFs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Mutect-GRCh37", full.names = TRUE)
  files1 <- files[1:2]
  list.of.vcfs1 <- ReadAndSplitMutectVCFs(files1)
  expect_equal(length(list.of.vcfs1$discarded.variants), 1)
  
  files2 <- files
  list.of.vcfs2 <- 
    expect_warning(ReadAndSplitMutectVCFs(files2,
                                          suppress.discarded.variants.warnings = FALSE))
  expect_equal(length(list.of.vcfs2$discarded.variants), 2)
    
  SBS.cats <- 
    expect_warning(VCFsToSBSCatalogs(list.of.vcfs2$SBS, ref.genome = "hg19", 
                                     region = "genome", 
                                     suppress.discarded.variants.warnings = FALSE))
  expect_false(is.null(SBS.cats$discarded.variants))
  SBS.cats1 <- 
    VCFsToSBSCatalogs(list.of.vcfs2$SBS, ref.genome = "hg19", 
                      region = "genome", return.annotated.vcfs = TRUE)
  expect_false(is.null(SBS.cats1$annotated.vcfs))
  
  DBS.cats <- 
    expect_warning(VCFsToDBSCatalogs(list.of.vcfs2$DBS, ref.genome = "hg19",
                                     region = "genome",
                                     suppress.discarded.variants.warnings = FALSE))
  expect_false(is.null(DBS.cats$discarded.variants))
  DBS.cats1 <- 
    VCFsToDBSCatalogs(list.of.vcfs2$DBS, ref.genome = "hg19",
                      region = "genome", return.annotated.vcfs = TRUE)
  expect_false(is.null(DBS.cats1$annotated.vcfs))
    
})

test_that("ReadAndSplitStrelkaSBSVCFs", {
  files <- list.files(path = "testdata/Strelka-SBS-GRCh37", full.names = TRUE)
  files1 <- files[1:2]
  list.of.vcfs1 <- ReadAndSplitStrelkaSBSVCFs(files1)
  expect_null(list.of.vcfs1$discarded.variants)
  
  files2 <- files
  list.of.vcfs2 <- 
    expect_warning(ReadAndSplitStrelkaSBSVCFs(files2,
                                              suppress.discarded.variants.warnings = FALSE))
  
  expect_false(is.null(list.of.vcfs2$discarded.variants))
})

test_that("ReadStrelkaIDVCFs", {
  files <- list.files(path = "testdata/Strelka-ID-GRCh37", full.names = TRUE)
  files1 <- files[1:2]
  list.of.vcfs1 <- ReadStrelkaIDVCFs(files1)
  expect_null(list.of.vcfs1$discarded.variants)
  
  files2 <- files
  list.of.vcfs2 <- ReadStrelkaIDVCFs(files2)
  
  expect_null(list.of.vcfs2$discarded.variants)
})

test_that("MutectVCFFilesToCatalog and VCFsToCatalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Mutect-GRCh37/", full.names = TRUE)
  catalogs <- MutectVCFFilesToCatalog(files, ref.genome = "hg19", 
                                      region = "genome")
  expect_false(is.null(catalogs$discarded.variants))
  expect_null(catalogs$annotated.vcfs)
  
  catalogs1 <- 
    MutectVCFFilesToCatalog(files, ref.genome = "hg19", 
                            region = "genome", return.annotated.vcfs = TRUE)
  expect_false(is.null(catalogs1$annotated.vcfs))
  
  catalogs2 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "mutect", region = "genome")
  expect_equal(catalogs, catalogs2)
  catalogs3 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "mutect", region = "genome",
                              return.annotated.vcfs = TRUE)
  expect_equal(catalogs1, catalogs3)
})

test_that("StrelkaSBSVCFFilesToCatalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Strelka-SBS-GRCh37", full.names = TRUE)
  catalogs <- StrelkaSBSVCFFilesToCatalog(files, ref.genome = "hg19", 
                                          region = "genome")
  expect_false(is.null(catalogs$discarded.variants))
  expect_null(catalogs$annotated.vcfs)
  
  catalogs1 <- 
    StrelkaSBSVCFFilesToCatalog(files, ref.genome = "hg19", 
                                region = "genome", return.annotated.vcfs = TRUE)
  expect_false(is.null(catalogs1$annotated.vcfs))
  
  catalogs2 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "strelka", region = "genome")
  catalogs2$catID <- NULL
  expect_equal(catalogs, catalogs2)
  
  catalogs3 <- VCFsToCatalogs(files, ref.genome = "hg19", variant.caller = "strelka",
                              region = "genome", return.annotated.vcfs = TRUE)
  catalogs3$catID <- catalogs3$annotated.vcfs$ID <- NULL
  expect_equal(catalogs1, catalogs3)
})

test_that("StrelkaIDVCFFilesToCatalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Strelka-ID-GRCh37", full.names = TRUE)
  catalogs <- StrelkaIDVCFFilesToCatalog(files, ref.genome = "hg19", 
                                         region = "genome")
  expect_false(is.null(catalogs$discarded.variants))
  expect_null(catalogs$annotated.vcfs)
  
  catalogs1 <- 
    StrelkaIDVCFFilesToCatalog(files, ref.genome = "hg19", 
                               region = "genome", return.annotated.vcfs = TRUE)
  expect_false(is.null(catalogs1$annotated.vcfs))
  
  catalogs2 <- VCFsToCatalogs(files, ref.genome = "hg19", 
                              variant.caller = "strelka", region = "genome")
  expect_equal(catalogs$catalog, catalogs2$catID)        
  expect_equal(catalogs$discarded.variants, catalogs2$discarded.variants)  
  
  catalogs3 <- VCFsToCatalogs(files, ref.genome = "hg19", variant.caller = "strelka",
                              region = "genome", return.annotated.vcfs = TRUE)
  expect_equal(catalogs1$catalog, catalogs3$catID)        
  expect_equal(catalogs1$discarded.variants, catalogs3$discarded.variants)  
  expect_equal(catalogs1$annotated.vcfs, catalogs3$annotated.vcfs$ID)  
  
})
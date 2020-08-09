context("Test functions dealing with discarded variants")

test_that("ReadAndSplitMutectVCFs", {
  files <- list.files(path = "testdata/Mutect-GRCh37/", full.names = TRUE)
  list.of.vcfs <- ReadAndSplitMutectVCFs(files)
  expect_null(list.of.vcfs$other.subs)
  expect_false(is.null(list.of.vcfs$multiple.alt))
  expect_null(list.of.vcfs$not.analyzed)
  
  file1 <- "testdata/Mutect.GRCh37.with.discarded.variants.vcf"
  files1 <- c(files, file1)
  list.of.vcfs1 <- ReadAndSplitMutectVCFs(files1)
  list.of.vcfs2 <- 
    expect_warning(ReadAndSplitMutectVCFs(files1,
                                          suppress.discarded.variants.warnings = FALSE))
    
  expect_false(is.null(list.of.vcfs1$multiple.alt))
  expect_false(is.null(list.of.vcfs1$other.subs))
  expect_false(is.null(list.of.vcfs1$not.analyzed))
})

test_that("MutectVCFFilesToCatalog", {
  files <- list.files(path = "testdata/Mutect-GRCh37/", full.names = TRUE)
  catalogs <- MutectVCFFilesToCatalog(files, ref.genome = "hg19", 
                                      region = "genome")
  expect_false(is.null(catalogs$discarded.variants))
  expect_null(catalogs$annotated.vcfs)
  
  file1 <- "testdata/Mutect.GRCh37.with.discarded.variants.vcf"
  files1 <- c(files, file1)
  catalogs1 <- MutectVCFFilesToCatalog(files, ref.genome = "hg19", 
                                       region = "genome")
  
  
})
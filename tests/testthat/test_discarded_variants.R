context("Test functions dealing with discarded variants and annotated vcfs")

test_that("ReadAndSplitMutectVCFs", {
  files <- list.files(path = "testdata/Mutect-GRCh37", full.names = TRUE)
  files1 <- files[1:2]
  list.of.vcfs1 <- ReadAndSplitMutectVCFs(files1)
  expect_null(list.of.vcfs1$other.subs)
  expect_false(is.null(list.of.vcfs1$multiple.alt))
  expect_null(list.of.vcfs1$not.analyzed)
  
  files2 <- files
  list.of.vcfs2 <- ReadAndSplitMutectVCFs(files2)
  list.of.vcfs3 <- 
    expect_warning(ReadAndSplitMutectVCFs(files2,
                                          suppress.discarded.variants.warnings = FALSE))
    
  expect_false(is.null(list.of.vcfs2$multiple.alt))
  expect_false(is.null(list.of.vcfs2$other.subs))
  expect_false(is.null(list.of.vcfs2$not.analyzed))
})

test_that("ReadAndSplitStrelkaSBSVCFs", {
  files <- list.files(path = "testdata/Strelka-SBS-GRCh37", full.names = TRUE)
  files1 <- files[1:2]
  list.of.vcfs1 <- ReadAndSplitStrelkaSBSVCFs(files1)
  expect_null(list.of.vcfs1$ThreePlus)
  expect_null(list.of.vcfs1$multiple.alt)
  expect_null(list.of.vcfs1$not.analyzed)
  
  files2 <- files
  list.of.vcfs2 <- ReadAndSplitStrelkaSBSVCFs(files2)
  list.of.vcfs3 <- 
    expect_warning(ReadAndSplitStrelkaSBSVCFs(files2,
                                              suppress.discarded.variants.warnings = FALSE))
  
  expect_null(list.of.vcfs2$ThreePlus)
  expect_false(is.null(list.of.vcfs2$multiple.alt))
  expect_false(is.null(list.of.vcfs2$not.analyzed))
})

test_that("MutectVCFFilesToCatalog", {
  files <- list.files(path = "testdata/Mutect-GRCh37/", full.names = TRUE)
  catalogs <- MutectVCFFilesToCatalog(files, ref.genome = "hg19", 
                                      region = "genome")
  expect_false(is.null(catalogs$discarded.variants))
  expect_null(catalogs$annotated.vcfs)
  
  catalogs1 <- 
    MutectVCFFilesToCatalog(files, ref.genome = "hg19", 
                            region = "genome", return.annotated.vcfs = TRUE)
  expect_false(is.null(catalogs1$annotated.vcfs))
})

test_that("StrelkaSBSVCFFilesToCatalog", {
  files <- list.files(path = "testdata/Strelka-SBS-GRCh37", full.names = TRUE)
  catalogs <- StrelkaSBSVCFFilesToCatalog(files, ref.genome = "hg19", 
                                          region = "genome")
  expect_false(is.null(catalogs$discarded.variants))
  expect_null(catalogs$annotated.vcfs)
  
  catalogs1 <- 
    StrelkaSBSVCFFilesToCatalog(files, ref.genome = "hg19", 
                                region = "genome", return.annotated.vcfs = TRUE)
  expect_false(is.null(catalogs1$annotated.vcfs))
})
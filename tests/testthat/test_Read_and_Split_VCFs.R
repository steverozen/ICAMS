context("Tests for internal functions reading and splitting vcfs")

test_that("Test ReadMutectVCFs and SplitListOfMutectVCFs", {
  files <- list.files(path = "testdata/Mutect-GRCh37", full.names = TRUE)
  vcfs <- ReadMutectVCFs(files)
  expect_equal(dim(vcfs[[1]]), c(1851,13))
  
  split.vcfs <- SplitListOfMutectVCFs(vcfs)
  expect_null(split.vcfs$discarded.variants$`Mutect.GRCh37.s2`)
})

test_that("Test ReadStrelkaSBSVCFs and SplitListOfStrelkaSBSVCFs", {
  files <- list.files(path = "testdata/Strelka-SBS-GRCh37/", full.names = TRUE)
  vcfs <- ReadStrelkaSBSVCFs(files)
  expect_equal(dim(vcfs[[1]]), c(798, 21))
  
  split.vcfs <- SplitListOfStrelkaSBSVCFs(vcfs)
  expect_null(split.vcfs$discarded.variants$`Mutect.GRCh37.s2`)
})
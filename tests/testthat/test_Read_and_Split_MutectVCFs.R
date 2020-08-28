context("Tests for internal functions reading and splitting mutect vcfs")

test_that("Test ReadMutectVCFs and SplitListOfMutectVCFs", {
  files <- list.files(path = "testdata/Mutect-GRCh37", full.names = TRUE)
  vcfs <- ReadMutectVCFs(files)
  expect_equal(dim(vcfs[[1]]), c(1851,13))
  
  split.vcfs <- SplitListOfMutectVCFs(vcfs)
  
})
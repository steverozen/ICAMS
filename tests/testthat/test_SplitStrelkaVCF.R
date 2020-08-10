context("Split Strelka SBS VCF")

test_that("StrelkaSBSVCFFilesToCatalog", {
  vcf <- ReadStrelkaSBSVCFs("testdata/Strelka-mult-alt.vcf")
  expect_equal(dim(vcf[[1]]$df), c(459, 13))
  split.vcfs <- expect_warning(SplitStrelkaSBSVCF(vcf[[1]]$df))
  expect_equal(dim(split.vcfs$multiple.alt), c(2, 13))
  expect_equal(dim(split.vcfs$SBS.vcf), c(447, 13))
})

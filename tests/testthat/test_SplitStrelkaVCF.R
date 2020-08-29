context("Split Strelka SBS VCF")

test_that("StrelkaSBSVCFFilesToCatalog", {
  vcf <- ReadStrelkaSBSVCFs("testdata/Strelka-mult-alt.vcf")
  expect_equal(dim(vcf[[1]]), c(459, 13))
  split.vcfs <- expect_warning(SplitStrelkaSBSVCF(vcf[[1]]))
  expect_equal(dim(split.vcfs$discarded.variants), c(2, 14))
  expect_equal(dim(split.vcfs$SBS.vcf), c(447, 13))
})

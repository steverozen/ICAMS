context("Split Strelka SBS VCF")

test_that("StrelkaSBSVCFFilesToCatalog", {
  vcf <- ReadStrelkaSBSVCFs("testdata/Strelka-mult-alt.vcf")
  expect_equal(dim(vcf[[1]]), c(459, 12))
  split.vcfs <- SplitStrelkaSBSVCF(vcf[[1]])
  expect_equal(dim(split.vcfs$multiple.alt), c(2, 12))
  expect_equal(dim(split.vcfs$SBS.vcf), c(447, 12))
})

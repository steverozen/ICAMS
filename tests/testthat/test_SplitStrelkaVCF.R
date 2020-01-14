context("Split Strelka SBS VCF")

test_that("StrelkaSBSVCFFilesToCatalog", {
  list <- ReadStrelkaSBSVCFs("testdata/Strelka-mult-alt.vcf")
  expect_equal(dim(list$`Strelka-mult-alt`[[1]]), c(459, 12))
  list1 <- SplitStrelkaSBSVCF(list$`Strelka-mult-alt`[[1]])
  expect_equal(dim(list1$multiple.alt), c(2, 12))
  expect_equal(dim(list1$SBS.vcf), c(447, 12))
})

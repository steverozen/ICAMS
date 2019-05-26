context("ReadStrelkaSBSVCF")

test_that('Read in VCF with column that has name "strand" or "VAF"', {
  expect_warning(ReadStrelkaSBSVCF("testdata/Strelka.colname.strand.vcf"))
  expect_warning(ReadStrelkaSBSVCF("testdata/Strelka.colname.VAF.vcf"))
})

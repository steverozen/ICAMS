context("ReadStrelkaSNSVCF")

test_that('Read in VCF with column that has name "strand" or "VAF"', {
  expect_warning(ReadStrelkaSNSVCF("testdata/Strelka.colname.strand.vcf"))
  expect_warning(ReadStrelkaSNSVCF("testdata/Strelka.colname.VAF.vcf"))
})
